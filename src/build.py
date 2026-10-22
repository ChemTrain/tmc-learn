"""Build TMC Learn.

Counts every complex in tmc_sets.py by both methods (tmc_count.py) and stops
if a count differs from the one typed in the set. Draws every complex with
RDKit: the metal is placed by hand, each ligand is drawn by RDKit from its
SMILES fragment and set on a site of a fixed geometry (octahedron, square
plane, piano stool, ...), and η-ligands get a line from the metal to the
centre of the bound atoms. The sets and drawings are injected into
template.html.

    python build.py                  # writes ../tmc_learn.html and dist/tmc_learn.html
    python build.py --sheets DIR     # also writes one PNG contact sheet per set (needs rsvg-convert)

../tmc_learn.html is the double-click version (a full HTML document).
dist/tmc_learn.html is the same page without <head>; it is what gets published
as the Claude artifact, whose host adds the document wrapper itself.
"""
import argparse
import base64
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point3D

from tmc_count import GROUP, count, fmt_os, parse_ligands
from tmc_ligands import LIGANDS, contribution
from tmc_sets import SETS

HERE = Path(__file__).resolve().parent
TEMPLATE = HERE / "template.html"
FRAGMENT = HERE / "dist" / "tmc_learn.html"
STANDALONE = HERE.parent / "tmc_learn.html"
PLACEHOLDER = "/*__TMC_DATA__*/null"
STANDALONE_HEAD = (
    "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, viewport-fit=cover\">\n"
)

LOGO = HERE.parent / "logo.svg"  # the ChemTrain logo; README.md shows the same file


def logo_parts():
    """The logo's drawing (without the outer <svg>) and its viewBox numbers."""
    svg = LOGO.read_text(encoding="utf-8")
    box = tuple(map(float, re.search(r'viewBox="([^"]+)"', svg).group(1).split()))
    inner = svg[svg.index(">", svg.index("<svg")) + 1:svg.rindex("</svg>")]
    return re.sub(r"<title>.*?</title>", "", inner), box


def add_logo(page):
    """Put the logo wherever the template has <!--__LOGO__-->; the drawing itself is included once."""
    if "<!--__LOGO__-->" not in page or "<!--__LOGO_DEFS__-->" not in page:
        sys.exit("template needs <!--__LOGO__--> and <!--__LOGO_DEFS__--> marks")
    inner, (x, y, w, h) = logo_parts()
    defs = (f'<svg width="0" height="0" style="position:absolute" aria-hidden="true">'
            f'<symbol id="ct-logo" viewBox="{x} {y} {w} {h}">{inner}</symbol></svg>')
    use = f'<svg class="logo" viewBox="0 0 {w} {h}" aria-hidden="true"><use href="#ct-logo"/></svg>'
    return page.replace("<!--__LOGO_DEFS__-->", defs).replace("<!--__LOGO__-->", use)


def favicon():
    """A square version of the logo for the browser tab of the stand-alone page."""
    inner, (x, y, w, h) = logo_parts()
    s = max(w, h)
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{x - (s - w) / 2:.2f} {y - (s - h) / 2:.2f} '
           f'{s:.2f} {s:.2f}">{inner}</svg>')
    return f'<link rel="icon" type="image/svg+xml" href="data:image/svg+xml;base64,{base64.b64encode(svg.encode()).decode()}">\n'


BOND = 26  # drawing units per bond; the page scales the SVG
# Atoms are drawn in placeholder colours, then swapped for CSS classes so the
# page can recolour them for light and dark themes. Metals get their own class.
COLOUR_CLASS = {1: "c", 2: "o", 3: "n", 4: "s", 5: "x", 6: "p", 7: "m"}
PALETTE = {-1: 1, 0: 1, 1: 1, 5: 1, 6: 1, 14: 1, 8: 2, 7: 3, 16: 4, 9: 5, 17: 5, 35: 5, 53: 5, 15: 6}
PALETTE.update({Chem.GetPeriodicTable().GetAtomicNumber(m): 7 for m in GROUP})
HEX_CLASS = {f"#0101{k:02X}": v for k, v in COLOUR_CLASS.items()}

rdDepictor.SetPreferCoordGen(True)
RDLogger.DisableLog("rdApp.*")

# ---------------------------------------------------------------- geometry
# A site is (angle in degrees, counter-clockwise from 3 o'clock, bond style[, length factor]):
# "" plain, "w" wedge (towards the viewer), "d" dash (away). Sites are filled
# in order; a chelate takes consecutive sites, so consecutive sites are cis.
SIGMA = 1.5  # metal–donor distance, in RDKit's bond units (1.5)

SITES = {
    "lin2": [(180, ""), (0, "")],
    "trig3": [(90, ""), (210, ""), (330, "")],
    "tet4": [(90, ""), (205, ""), (335, "w"), (270, "d")],
    "sqpl4": [(180, ""), (90, ""), (0, ""), (270, "")],
    "tbp5": [(90, ""), (0, ""), (215, "w"), (145, "d"), (270, "")],
    # square pyramid: apex, then the base front-left, back-left, back-right, front-right
    # (so sites 2/4 and 3/5 are trans, and a chelate on apex + site 2 spans a wide angle)
    "sqpy5": [(90, ""), (215, "w"), (160, "d"), (20, "d"), (325, "w")],
    # octahedron: up, front-left, front-right, back-right, back-left, down
    # trans pairs: 1/6, 2/4, 3/5; consecutive sites are always cis
    "oct6": [(90, ""), (210, "w"), (330, "w"), (30, "d"), (150, "d"), (270, "")],
    # the same octahedron for two trans chelates: up, the left pair, the right pair, down
    "oct6b": [(90, ""), (135, "d"), (225, "w"), (315, "w"), (45, "d"), (270, "")],
    # three chelates drawn flat, as a propeller (consecutive pairs are the chelates)
    "hex6": [(30, ""), (90, ""), (150, ""), (210, ""), (270, ""), (330, "")],
    # two η²,η²-dienes (Ni(cod)₂), one on each side
    "tet22": [(150, ""), (210, ""), (30, ""), (330, "")],
    # piano stools and sandwiches: the ring first
    "ps2": [(90, ""), (225, ""), (315, "")],
    "ps3": [(90, ""), (215, ""), (270, "w"), (325, "")],
    "ps4": [(90, ""), (195, ""), (240, "w"), (300, "w"), (345, "")],
    "sand": [(90, ""), (270, "")],
    # a flat tetradentate ligand (salen, porphyrin), drawn face-on, and one axial ligand below
    "pl4": [(225, ""), (135, ""), (45, ""), (315, "")],
    # one ligand on a generic metal M, for the Ligands set
    "L1": [(90, "")],
    "L2": [(135, ""), (45, "")],
    "L3": [(180, ""), (90, ""), (0, "")],
    "pl4ax": [(225, ""), (135, ""), (45, ""), (315, ""), (270, "")],
    "pl4axl": [(225, ""), (135, ""), (45, ""), (315, ""), (270, "", 1.6)],
    "bent1": [(135, ""), (225, ""), (0, "")],
    "bent2": [(135, ""), (225, ""), (35, ""), (325, "")],
    "bent3": [(135, ""), (225, ""), (45, ""), (0, ""), (315, "")],
}


def unit(deg):
    a = math.radians(deg)
    return math.cos(a), math.sin(a)


def layout(geom):
    """Metal positions, metal–metal bonds, bridge sites and per-metal terminal sites.

    Returns (metals [(x, y)], mm [(i, j)], bridges [(x, y)], terminal [(metal, angle, style)]).
    """
    if geom in SITES:
        return [(0.0, 0.0)], [], [], [(0, *site) for site in SITES[geom]]
    a = SIGMA * math.cos(math.radians(45))
    if geom in ("d-sqpl", "d-sqpl1", "d-ps"):
        # two metals sharing an edge of two μ-ligands, no M–M bond
        metals, mm, bridges = [(-a, 0.0), (a, 0.0)], [], [(0.0, a), (0.0, -a)]
        term = {"d-sqpl": [(0, 135, ""), (0, 225, ""), (1, 45, ""), (1, 315, "")],
                "d-sqpl1": [(0, 180, ""), (1, 0, "")],
                "d-ps": [(0, 160, ""), (0, 235, "w"), (1, 20, ""), (1, 305, "w")]}[geom]
        return metals, mm, bridges, term
    if geom == "mu1":
        # one μ-ligand between two generic metals, for the Ligands set
        return [(-a, 0.0), (a, 0.0)], [], [(0.0, a)], []
    if geom == "d-mm5":
        d = 0.95
        metals = [(-d, 0.0), (d, 0.0)]
        term = [(0, 90, ""), (0, 220, "w"), (0, 180, ""), (0, 140, "d"), (0, 270, ""),
                (1, 90, ""), (1, 320, "w"), (1, 0, ""), (1, 40, "d"), (1, 270, "")]
        return metals, [(0, 1)], [], term
    if geom in ("d-co2", "d-fp"):
        d = 0.95
        metals, bridges = [(-d, 0.0), (d, 0.0)], [(0.0, 1.1), (0.0, -1.1)]
        term = {"d-co2": [(0, 180, ""), (0, 220, "w"), (0, 140, "d"), (1, 0, ""), (1, 320, "w"), (1, 40, "d")],
                "d-fp": [(0, 180, ""), (0, 235, "w"), (1, 0, ""), (1, 55, "w")]}[geom]
        return metals, [(0, 1)], bridges, term
    if geom == "tri":
        r = 1.1
        metals = [(r * math.cos(math.radians(t)), r * math.sin(math.radians(t))) for t in (90, 210, 330)]
        term = []
        for i, t in enumerate((90, 210, 330)):
            term += [(i, t + 55, ""), (i, t + 18, "d"), (i, t - 18, "w"), (i, t - 55, "")]
        return metals, [(0, 1), (1, 2), (0, 2)], [], term
    raise ValueError(f"unknown geometry {geom!r}")


# ---------------------------------------------------------------- ligand fragments
# geometries whose chelates sit further out than SIGMA, to leave room between them
REACH = {"hex6": 0.9, "oct6b": 0.45}


def coords(mol, coordgen=True):
    """2D coordinates; CoordGen gives NaN for molecules without heavy atoms (H₂)."""
    rdDepictor.SetPreferCoordGen(coordgen)
    rdDepictor.Compute2DCoords(mol)
    rdDepictor.SetPreferCoordGen(True)
    xy = mol.GetConformer().GetPositions()
    if any(math.isnan(v) for v in xy.ravel()):
        rdDepictor.SetPreferCoordGen(False)
        rdDepictor.Compute2DCoords(mol)
        rdDepictor.SetPreferCoordGen(True)
        xy = mol.GetConformer().GetPositions()
    return [(float(p[0]), float(p[1])) for p in xy]


def fragment(key):
    """Parse a ligand's SMILES and lay it out.

    Returns (mol, xy [(x, y)], points [[atom indices]], labels {atom index: (right, left)}, mpos),
    where mpos is where the metal sits in xy, or None when the ligand binds at one point
    and is simply turned to face the metal.
    """
    lig = LIGANDS[key]
    params = Chem.SmilesParserParams()
    params.removeHs = False
    params.sanitize = False
    mol = Chem.MolFromSmiles(lig["smiles"], params)
    if mol is None:
        raise ValueError(f"RDKit cannot parse the {key} SMILES")
    mol.UpdatePropertyCache(strict=False)
    by_map = {a.GetAtomMapNum(): a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum()}
    donors = sorted(k for k in by_map if k < 50)
    points = [[by_map[k] for k in p] for p in lig.get("points", [(k,) for k in donors])]
    labels = {by_map[k]: v for k, v in lig.get("labels", {}).items()}
    for a in mol.GetAtoms():
        a.SetAtomMapNum(0)
    mpos = None
    if "xy" in lig:
        xy = [tuple(map(float, p)) for p in lig["xy"]]
        if len(xy) != mol.GetNumAtoms():
            raise ValueError(f"{key}: {len(xy)} coordinates for {mol.GetNumAtoms()} atoms")
        mpos = tuple(map(float, lig["mpos"]))
    elif len(points) > 1 and all(len(p) == 1 for p in points):
        # a chelate: lay it out bound to a stand-in metal, so the chelate ring comes out right
        rw = Chem.RWMol(mol)
        dummy = rw.AddAtom(Chem.Atom(0))
        for p, order in zip(points, lig.get("order", [1] * len(points))):
            rw.AddBond(dummy, p[0], BOND_TYPE[order])
        bound = rw.GetMol()
        bound.UpdatePropertyCache(strict=False)
        xy = coords(bound, lig.get("depictor") != "rdkit")
        mpos = xy.pop(dummy)
    elif mol.GetNumAtoms() == 1:
        xy = [(0.0, 0.0)]
    else:
        xy = coords(mol)
    return mol, xy, points, labels, mpos


def mean(pts):
    return sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts)


def transform(xy, rot, mirror, shift):
    c, s = math.cos(rot), math.sin(rot)
    out = []
    for x, y in xy:
        if mirror:
            y = -y
        out.append((x * c - y * s + shift[0], x * s + y * c + shift[1]))
    return out


def place_one(xy, point, centre, direction):
    """Put a single-point ligand on a site: the point at `centre`, the rest pointing along `direction`."""
    p = mean([xy[i] for i in point])
    rel = [(x - p[0], y - p[1]) for x, y in xy]
    want = math.atan2(direction[1], direction[0])
    if len(point) == 2:
        # η²: the bound bond lies across the M–L direction, the rest of the ligand points away
        a, b = rel[point[0]], rel[point[1]]
        rot = want + math.pi / 2 - math.atan2(b[1] - a[1], b[0] - a[0])
        out = transform(rel, rot, False, centre)
        com = mean(out)
        if (com[0] - centre[0]) * direction[0] + (com[1] - centre[1]) * direction[1] < -1e-6:
            out = transform(rel, rot, True, centre)
            # mirroring flips the bound bond; turn it back across the M–L direction
            a2, b2 = out[point[0]], out[point[1]]
            fix = want + math.pi / 2 - math.atan2(b2[1] - a2[1], b2[0] - a2[0])
            out = transform([(x - centre[0], y - centre[1]) for x, y in out], fix, False, centre)
        return out
    com = mean(rel)
    if math.hypot(*com) > 0.05:
        rot = want - math.atan2(com[1], com[0])
    elif len(point) > 2:
        # a bare ring: point one ring atom away from the metal
        r0 = rel[point[0]]
        rot = want - math.atan2(r0[1], r0[0])
    else:
        rot = 0.0
    return transform(rel, rot, False, centre)


def place_pivot(xy, points, mpos, metal, angles):
    """Turn a ligand laid out around a stand-in metal (mpos) so its points face the sites at `angles`.

    Tries the ligand and its mirror image; keeps whichever points match the sites best.
    """
    vecs = []
    for p in points:
        c = mean([xy[i] for i in p])
        d = math.hypot(c[0] - mpos[0], c[1] - mpos[1]) or 1
        vecs.append(((c[0] - mpos[0]) / d, (c[1] - mpos[1]) / d))
    rel = [(x - mpos[0], y - mpos[1]) for x, y in xy]
    best = None
    for mirror in (False, True):
        a = [(x, -y if mirror else y) for x, y in vecs]
        b = [unit(t) for t in angles]
        num = sum(p[0] * q[1] - p[1] * q[0] for p, q in zip(a, b))
        den = sum(p[0] * q[0] + p[1] * q[1] for p, q in zip(a, b))
        rot = math.atan2(num, den)
        c, s_ = math.cos(rot), math.sin(rot)
        fit = sum((p[0] * c - p[1] * s_) * q[0] + (p[0] * s_ + p[1] * c) * q[1] for p, q in zip(a, b))
        if best is None or fit > best[0] + 1e-6:
            best = (fit, rot, mirror)
    return transform(rel, best[1], best[2], metal)


def hapto_distance(xy, point):
    """How far the centre of an η-bound group sits from the metal."""
    if len(point) <= 2:
        return SIGMA
    c = mean([xy[i] for i in point])
    return max(math.hypot(xy[i][0] - c[0], xy[i][1] - c[1]) for i in point) + 1.0


BOND_TYPE = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
BOND_DIR = {"w": Chem.BondDir.BEGINWEDGE, "d": Chem.BondDir.BEGINDASH}


def assemble(metal, ligands, geom):
    """Build the drawing molecule of a complex.

    Returns (mol, metal atom indices, haptic lines [(metal idx, [atom idx])],
    anchors [(metal idx, [atom idx]) per ligand instance]).
    """
    metals, mm, bridges, term = layout(geom)
    sym = re.match(r"[A-Z][a-z]?", metal).group(0)
    if len(metals) != (int(metal[len(sym):] or 1)):
        raise ValueError(f"geometry {geom} has {len(metals)} metals, the formula {metal}")
    rw = Chem.RWMol()
    pos = {}
    midx = []
    for x, y in metals:
        # "M" is a generic metal: drawn as a metal atom (so it gets the metal colour), labelled M
        i = rw.AddAtom(Chem.Atom("Ti" if sym == "M" else sym))
        if sym == "M":
            rw.GetAtomWithIdx(i).SetProp("atomLabel", "M")
        pos[i] = (x, y)
        midx.append(i)
    for i, j in mm:
        rw.AddBond(midx[i], midx[j], Chem.BondType.SINGLE)

    centre = mean(metals)
    hapto, anchors, label_of = [], [], {}
    bi = ti = 0
    for key in ligands:
        lig = LIGANDS[key]
        if not lig.get("smiles"):
            raise ValueError(f"{key} has no drawing")
        mol, xy, points, labels, mpos = fragment(key)
        orders = lig.get("order", [1] * len(points))
        if lig.get("bridging"):
            if bi >= len(bridges):
                raise ValueError(f"geometry {geom} has only {len(bridges)} bridge sites")
            b = bridges[bi]
            bi += 1
            d = (b[0] - centre[0], b[1] - centre[1])
            n = math.hypot(*d) or 1
            out = place_one(xy, points[0], b, (d[0] / n, d[1] / n))
            bonded = [(m, points[0], 1, "") for m in range(len(metals))]
        else:
            k = len(points)
            if ti + k > len(term):
                raise ValueError(f"geometry {geom} has only {len(term)} sites")
            sites = term[ti:ti + k]
            ti += k
            if len({site[0] for site in sites}) != 1:
                raise ValueError(f"{key} would span two metals")
            m = sites[0][0]
            mx, my = metals[m]
            if mpos is not None:
                out = place_pivot(xy, points, mpos, (mx, my), [site[1] for site in sites])
                if geom in REACH:
                    u = mean([unit(site[1]) for site in sites])
                    n = math.hypot(*u) or 1
                    out = [(x + u[0] / n * REACH[geom], y + u[1] / n * REACH[geom]) for x, y in out]
            elif k == 1:
                u = unit(sites[0][1])
                dist = hapto_distance(xy, points[0]) * (sites[0][3] if len(sites[0]) > 3 else 1)
                out = place_one(xy, points[0], (mx + u[0] * dist, my + u[1] * dist), u)
            else:
                raise ValueError(f"{key} binds at {k} points but has no layout around a metal")
            bonded = [(m, p, o, site[2]) for site, p, o in zip(sites, points, orders)]

        off = rw.GetNumAtoms()
        for a in mol.GetAtoms():
            rw.AddAtom(a)
        for b in mol.GetBonds():
            nb = rw.AddBond(b.GetBeginAtomIdx() + off, b.GetEndAtomIdx() + off, b.GetBondType()) - 1
            if b.GetIsAromatic():
                rw.GetBondWithIdx(nb).SetIsAromatic(True)
        for i, p in enumerate(out):
            pos[off + i] = p
        for i, (right, left) in labels.items():
            label_of[off + i] = (right, left)
        bound = []
        for m, p, o, style in bonded:
            atoms = [off + i for i in p]
            bound += atoms
            if len(atoms) == 1:
                rw.AddBond(midx[m], atoms[0], BOND_TYPE[o])
                bond = rw.GetBondBetweenAtoms(midx[m], atoms[0])
                if style and o == 1:
                    bond.SetBondDir(BOND_DIR[style])
            else:
                hapto.append((midx[m], atoms))
        # (the metals it binds, the atoms that bind, all its atoms): draw() puts the badge past its far end
        anchors.append((sorted({midx[m] for m, _, _, _ in bonded}), bound, list(range(off, off + mol.GetNumAtoms()))))
    if ti != len(term) or bi != len(bridges):
        raise ValueError(f"{len(ligands)} ligands leave sites of {geom} empty ({ti}/{len(term)}, {bi}/{len(bridges)})")

    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)
    # labels read away from the metal: "PPh3" on the right, "Ph3P" on the left
    for i, (right, left) in label_of.items():
        x = pos[i][0]
        nearest = min(midx, key=lambda m: math.dist(pos[m], pos[i]))
        mol.GetAtomWithIdx(i).SetProp("atomLabel", left if x < pos[nearest][0] - 0.3 else right)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, (x, y) in pos.items():
        conf.SetAtomPosition(i, Point3D(x, y, 0))
    mol.RemoveAllConformers()
    mol.AddConformer(conf)
    return mol, midx, hapto, anchors


def clashes(mol, min_dist=0.55):
    """Pairs of atoms that are not bonded but sit closer than min_dist (bond length 1.5)."""
    xy = mol.GetConformer().GetPositions()
    out = []
    for i in range(mol.GetNumAtoms()):
        for j in range(i + 1, mol.GetNumAtoms()):
            if mol.GetBondBetweenAtoms(i, j) is None and math.dist(xy[i][:2], xy[j][:2]) < min_dist:
                out.append((i, j))
    return out


# ---------------------------------------------------------------- drawing
def set_options(o, bond=BOND):
    o.clearBackground = False
    o.fixedBondLength = bond
    o.padding = 0.12
    o.baseFontSize = 0.75
    o.addStereoAnnotation = False
    o.prepareMolsBeforeDrawing = False
    o.flagCloseContactsDist = -1  # build.py reports clashes itself (see clashes())


def minify_path(d):
    d = re.sub(r"\s+", " ", d).strip()
    d = re.sub(r"([MLQCZ]) ", r"\1", d)
    d = re.sub(r" ([MLQCZ])", r"\1", d)
    d = d.replace(", ", ",")
    return re.sub(r"(\d)\.0(?=[ ,MLQCZ]|$)", r"\1", d)


def svg_body(svg, highlight=()):
    """RDKit's SVG paths with their colours turned into classes.

    b: a bond line, t: text or a filled wedge, h: a thin line (the hashes of a dashed bond);
    then the colour class (c, o, n, s, x, p, m).
    """
    head = re.search(r"viewBox='0 0 ([\d.]+) ([\d.]+)'>", svg)
    width, height = float(head.group(1)), float(head.group(2))
    body = svg[head.end():svg.rindex("</svg>")].replace("<!-- END OF HEADER -->", "")
    elements = re.findall(r"<[a-z]+\b[^>]*/>", body, flags=re.S)
    leftover = re.sub(r"<[a-z]+\b[^>]*/>", "", body, flags=re.S).strip()
    if leftover:
        raise ValueError(f"unexpected SVG content: {leftover[:80]}")
    out = []
    for el in elements:
        if not el.startswith("<path"):
            raise ValueError(f"unexpected SVG element: {el[:60]}")
        classes = (re.search(r"class='([^']*)'", el) or [None, ""])[1].split()
        path = minify_path(re.search(r"d='([^']*)'", el, flags=re.S).group(1))
        style = re.search(r"style='([^']*)'", el)
        stroke = style and re.search(r"stroke:(#[0-9A-Fa-f]{6})", style.group(1))
        sfill = style and re.search(r"fill:(#[0-9A-Fa-f]{6})", style.group(1))
        fill = re.search(r"fill='(#[0-9A-Fa-f]{6})'", el)
        if stroke and "fill:none" in style.group(1):
            width_px = float(re.search(r"stroke-width:([\d.]+)px", style.group(1)).group(1))
            kind, colour = ("h" if width_px < 1.5 else "b"), stroke.group(1)
        elif fill and not style:
            kind, colour = "t", fill.group(1)
        elif sfill:
            kind, colour = "t", sfill.group(1)
        else:
            raise ValueError(f"unexpected path style: {el[:120]}")
        cls = [kind, HEX_CLASS[colour.upper()]]
        if any(f"atom-{h}" in classes for h in highlight):
            cls.append("hh")
        out.append(f'<path class="{" ".join(cls)}" d="{path}"/>')
    return round(width, 1), round(height, 1), out


LABEL_GAP = 8.5  # drawing units kept clear around a metal symbol by a haptic line
BADGE_GAP = 9.5  # from a ligand's far end to the centre of its Count it badge
BADGE_CHAR = 4.6  # rough width of one label character, in drawing units


def draw(mol, midx, hapto, anchors):
    """Return (width, height, metal positions, ligand anchors, svg body) for the page.

    anchors hold [x, y, dx, dy] for each ligand: a point just past the ligand's
    far end, seen from the metal, where the page puts a letter badge, and the
    unit vector from the metal.
    """
    drawer = rdMolDraw2D.MolDraw2DSVG(-1, -1)
    o = drawer.drawOptions()
    set_options(o)
    o.setAtomPalette({k: (1 / 255, 1 / 255, v / 255) for k, v in PALETTE.items()})
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    width, height, paths = svg_body(drawer.GetDrawingText())

    def at(i):
        p = drawer.GetDrawCoords(i)
        return p.x, p.y

    for m, atoms in hapto:
        mx, my = at(m)
        cx, cy = mean([at(i) for i in atoms])
        d = math.hypot(cx - mx, cy - my) or 1
        sx, sy = mx + (cx - mx) / d * LABEL_GAP, my + (cy - my) / d * LABEL_GAP
        # stop short of a ring centre, so the line meets the ring rather than crossing it
        stop = 0 if len(atoms) == 2 else 4
        ex, ey = cx - (cx - mx) / d * stop, cy - (cy - my) / d * stop
        paths.insert(0, f'<path class="b c" d="M{sx:.1f},{sy:.1f}L{ex:.1f},{ey:.1f}"/>')

    metals = [[round(v, 1) for v in at(m)] for m in midx]
    out = []
    for ms, bound, atoms in anchors:
        ox, oy = mean([at(m) for m in ms])
        cx, cy = mean([at(i) for i in bound])
        d = math.hypot(cx - ox, cy - oy) or 1
        ux, uy = (cx - ox) / d, (cy - oy) / d
        reach = 0
        for i in atoms:
            x, y = at(i)
            ext = 0
            atom = mol.GetAtomWithIdx(i)
            if atom.HasProp("atomLabel"):  # a text label (PPh₃, Mes) sticks out past its atom
                text = re.sub(r"<[^>]+>", "", atom.GetProp("atomLabel"))
                ext = len(text) * BADGE_CHAR * abs(ux) + 5 * abs(uy)
            elif atom.GetAtomicNum() not in (0, 6) or not atom.GetDegree():
                # an element symbol, with its H's (Cl, NH₃)
                chars = len(atom.GetSymbol()) + (2 if atom.GetTotalNumHs() else 0)
                ext = chars * BADGE_CHAR * abs(ux) + 5 * abs(uy)
            reach = max(reach, (x - ox) * ux + (y - oy) * uy + ext)
        reach += BADGE_GAP
        out.append([round(ox + ux * reach, 1), round(oy + uy * reach, 1), round(ux, 2), round(uy, 2)])
    return width, height, metals, out, "".join(paths)


# ---------------------------------------------------------------- data
def ligand_table():
    """Everything the page needs to show a ligand's contribution, keyed like LIGANDS."""
    out = {}
    for key, lig in LIGANDS.items():
        cov, ion_e, ion_q, nx = contribution(key)
        out[key] = {"name": lig["name"], "show": lig["show"], "cls": lig["cls"], "cov": cov,
                    "ion": ion_e, "q": ion_q, "ionName": lig.get("ion", ""), "mu": bool(lig.get("bridging"))}
    return out


# the real geometry each drawing layout stands for, for the Geometry questions; layouts with
# rings (piano stools, sandwiches), several metals or no drawing have none
SHAPE_OF_GEOM = {"lin2": "lin", "trig3": "trig", "tet4": "tet", "tet22": "tet", "sqpl4": "sqpl", "pl4": "sqpl",
                 "tbp5": "tbp", "sqpy5": "sqpy", "pl4ax": "sqpy", "pl4axl": "sqpy",
                 "oct6": "oct", "oct6b": "oct", "hex6": "oct"}
SITES_OF_SHAPE = {"lin": 2, "trig": 3, "tet": 4, "sqpl": 4, "tbp": 5, "sqpy": 5, "oct": 6, "tprism": 6}


def complex_card(c):
    ligs = parse_ligands(c["ligands"])
    r = count(c["metal"], c["charge"], ligs, c["mm"])
    problems = []
    if r["ve"] != c["ve"]:
        problems.append(f"electron count typed {c['ve']}, computed {r['ve']}")
    if r["os"] != c["os"]:
        problems.append(f"oxidation state typed {fmt_os(c['os'])}, computed {fmt_os(r['os'])}")
    if problems:
        raise ValueError("; ".join(problems))
    shape = c["shape"] if c["shape"] is not None else SHAPE_OF_GEOM.get(c["geom"], "")
    if shape and shape not in SITES_OF_SHAPE:
        raise ValueError(f"unknown shape {shape!r}")
    grouped, seq = [], []
    for k in ligs:  # [key, how many], in order of first appearance; seq: the group of each ligand, in drawing order
        for i, g in enumerate(grouped):
            if g[0] == k:
                g[1] += 1
                seq.append(i)
                break
        else:
            seq.append(len(grouped))
            grouped.append([k, 1])
    return {"term": c["name"], "formula": c["formula"], "metal": r["metal"], "n": r["n"], "charge": c["charge"],
            "mm": c["mm"], "ve": r["ve"], "os": r["os"], "dn": r["dn"], "group": r["group"], "lig": grouped, "seq": seq,
            "noOs": c["no_os"], "note": c["note"], "fam": c["fam"], "shape": shape}


def build_data():
    mols, index, problems, drawn = [], {}, [], {}
    warnings = []
    sets = []
    for s in SETS:
        cards = []
        for c in s["cards"]:
            if s["kind"] == "ligand":
                if c not in LIGANDS:
                    problems.append(f"{s['id']}: unknown ligand {c}")
                    continue
                card = {"term": LIGANDS[c]["name"], "formula": LIGANDS[c]["show"], "key": c, "mol": None}
                if LIGANDS[c].get("smiles"):
                    try:
                        spec = ("ligand", c)
                        if spec not in index:
                            if LIGANDS[c].get("bridging"):
                                geom, metal = "mu1", "M2"
                            else:
                                k = len(fragment(c)[2])
                                geom, metal = {1: "L1", 2: "L2", 3: "L3", 4: "pl4"}[k], "M"
                            mol, midx, hapto, anchors = assemble(metal, [c], geom)
                            width, height, metals, anc, body = draw(mol, midx, hapto, anchors)
                            index[spec] = len(mols)
                            mols.append([width, height, metals, anc, body])
                        card["mol"] = index[spec]
                    except Exception as exc:
                        problems.append(f"{s['id']}: {c}: {exc}")
                        continue
                cards.append(card)
                continue
            try:
                card = complex_card(c)
                card["mol"] = None
                if c["geom"]:
                    spec = (c["metal"], c["ligands"], c["geom"])
                    if spec not in index:
                        mol, midx, hapto, anchors = assemble(c["metal"], parse_ligands(c["ligands"]), c["geom"])
                        if clashes(mol):
                            warnings.append(f"{s['id']}: {c['name']}: {len(clashes(mol))} atom pairs overlap")
                        width, height, metals, anc, body = draw(mol, midx, hapto, anchors)
                        index[spec] = len(mols)
                        mols.append([width, height, metals, anc, body])
                        drawn[index[spec]] = (c["name"], width, height, metals, body)
                    card["mol"] = index[spec]
            except Exception as exc:  # report every bad card, not just the first
                problems.append(f"{s['id']}: {c['name']}: {exc}")
                continue
            cards.append(card)
        sets.append({k: s[k] for k in ("id", "name", "kind", "desc", "source")} | {"cards": cards})
    return {"sets": sets, "mols": mols, "ligands": ligand_table()}, problems, warnings, drawn


# ---------------------------------------------------------------- contact sheets
SHEET_CSS = """
.b{fill:none;stroke-width:1.35;stroke-linecap:round;stroke-linejoin:round}.h{fill:none;stroke-width:.8}
.b.c,.h.c{stroke:#1B2731}.b.o,.h.o{stroke:#C8261A}.b.n,.h.n{stroke:#2045C4}.b.s,.h.s{stroke:#9A7500}
.b.x,.h.x{stroke:#2A7A2F}.b.p,.h.p{stroke:#BF560B}.b.m,.h.m{stroke:#0D5C78}
.t.c{fill:#1B2731}.t.o{fill:#C8261A}.t.n{fill:#2045C4}.t.s{fill:#9A7500}.t.x{fill:#2A7A2F}.t.p{fill:#BF560B}.t.m{fill:#0D5C78}
.halo{fill:#F9CFE1}.anc{fill:#C0175A}.cap{font:11px sans-serif;fill:#333}
"""


def write_sheets(outdir, data, drawn):
    exe = shutil.which("rsvg-convert") or "/home/franci_sala/anaconda3/envs/aixchem-io/bin/rsvg-convert"
    outdir.mkdir(parents=True, exist_ok=True)
    cell_w, cell_h, cols = 230, 200, 6
    for s in data["sets"]:
        cells = [(c, data["mols"][c["mol"]]) for c in s["cards"] if c.get("mol") is not None]
        if not cells:
            continue
        rows = math.ceil(len(cells) / cols)
        parts = []
        for k, (c, (w, h, metals, anc, body)) in enumerate(cells):
            x0, y0 = (k % cols) * cell_w, (k // cols) * cell_h
            sc = min(1.6, (cell_w - 10) / w, (cell_h - 34) / h)
            g = "".join(f'<circle class="halo" cx="{x}" cy="{y}" r="9"/>' for x, y in metals)
            g += "".join(f'<circle class="anc" cx="{a[0]}" cy="{a[1]}" r="1.6"/>' for a in anc)
            ox, oy = x0 + (cell_w - w * sc) / 2, y0 + 4 + (cell_h - 34 - h * sc) / 2
            parts.append(f'<g transform="translate({ox:.1f},{oy:.1f}) scale({sc:.3f})">{g}{body}</g>')
            cap = f"{c['term'][:34]} · {c['ve']} e" if "ve" in c else c["term"][:40]
            parts.append(f'<text class="cap" x="{x0 + cell_w / 2}" y="{y0 + cell_h - 12}" text-anchor="middle">{cap}</text>')
        svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{cols * cell_w}" height="{rows * cell_h}">'
               f'<style>{SHEET_CSS}</style><rect width="100%" height="100%" fill="#fff"/>{"".join(parts)}</svg>')
        (outdir / f"{s['id']}.svg").write_text(svg, encoding="utf-8")
        subprocess.run([exe, str(outdir / f"{s['id']}.svg"), "-o", str(outdir / f"{s['id']}.png")], check=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sheets", type=Path, help="write PNG contact sheets of every set here")
    args = ap.parse_args()

    data, problems, warnings, drawn = build_data()
    if warnings:
        print("Drawings to check:", *warnings, sep="\n  ", file=sys.stderr)
    if problems:
        print("Problems:", *problems, sep="\n  ", file=sys.stderr)
        sys.exit(1)
    if args.sheets:
        write_sheets(args.sheets, data, drawn)

    n_cards = sum(len(s["cards"]) for s in data["sets"])
    n_named = sum(1 for s in data["sets"] for c in s["cards"] if c.get("fam") and c.get("mol") is not None)
    blob = json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
    if not TEMPLATE.exists():
        print(f"{len(data['sets'])} sets, {n_cards} cards, {len(data['mols'])} structures, "
              f"{n_named} named; data {len(blob.encode()) / 1024:.0f} KB (no template yet)")
        return
    template = TEMPLATE.read_text(encoding="utf-8")
    if PLACEHOLDER not in template:
        sys.exit(f"{TEMPLATE.name} has no {PLACEHOLDER} placeholder")
    page = add_logo(template.replace(PLACEHOLDER, blob))
    FRAGMENT.parent.mkdir(exist_ok=True)
    FRAGMENT.write_text(page, encoding="utf-8")
    STANDALONE.write_text(STANDALONE_HEAD + favicon() + page, encoding="utf-8")
    print(f"{len(data['sets'])} sets, {n_cards} cards, {len(data['mols'])} structures, "
          f"{n_named} named, page {len(page.encode()) / 1024:.0f} KB")
    print(f"wrote {STANDALONE}\nwrote {FRAGMENT}")


if __name__ == "__main__":
    main()
