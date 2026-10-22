"""Electron counting for TMC Learn.

count() works out the oxidation state, d-electron count and valence electron
count of a complex from its metal, overall charge, ligand list and number of
metal–metal bonds, by both methods:

    covalent   n·group + Σ(2L + X) − charge + 2·(M–M bonds), per metal ÷ n
    ionic      dⁿ + Σ ionic electrons ÷ n + 2·(M–M bonds) ÷ n,
               with oxidation state = (charge − Σ ionic ligand charges) ÷ n

and refuses a complex where the two disagree, or where a per-metal value is
not a whole number (a cluster that the simple per-metal average cannot count).
A μ-ligand is counted once for the whole molecule; each M–M bond gives one
electron to each of its two metals and does not change the oxidation state.
"""
import re

from tmc_ligands import LIGANDS, contribution

GROUP = {
    "Sc": 3, "Ti": 4, "V": 5, "Cr": 6, "Mn": 7, "Fe": 8, "Co": 9, "Ni": 10, "Cu": 11, "Zn": 12,
    "Y": 3, "Zr": 4, "Nb": 5, "Mo": 6, "Tc": 7, "Ru": 8, "Rh": 9, "Pd": 10, "Ag": 11, "Cd": 12,
    "Hf": 4, "Ta": 5, "W": 6, "Re": 7, "Os": 8, "Ir": 9, "Pt": 10, "Au": 11, "Hg": 12,
}


def parse_metal(text):
    """'Fe' -> ('Fe', 1); 'Ru3' -> ('Ru', 3)."""
    m = re.fullmatch(r"([A-Z][a-z]?)(\d*)", text)
    if not m or m.group(1) not in GROUP:
        raise ValueError(f"unknown metal {text!r}")
    return m.group(1), int(m.group(2) or 1)


def parse_ligands(text):
    """'CO*4 PPh3' -> ['CO', 'CO', 'CO', 'CO', 'PPh3'] (order kept: it is the drawing order)."""
    out = []
    for tok in text.split():
        key, _, times = tok.partition("*")
        if key not in LIGANDS:
            raise ValueError(f"unknown ligand {key!r}")
        out += [key] * int(times or 1)
    return out


def count(metal, charge, ligands, mm=0):
    sym, n = parse_metal(metal)
    group = GROUP[sym]
    if n == 1 and mm:
        raise ValueError("a metal–metal bond needs two metals")
    if n == 1 and any(LIGANDS[k].get("bridging") for k in ligands):
        raise ValueError("a bridging ligand needs two metals")
    cov = ion_e = ion_q = 0
    for k in ligands:
        c, e, q, _ = contribution(k)
        cov, ion_e, ion_q = cov + c, ion_e + e, ion_q + q

    ve_total = n * group + cov - charge + 2 * mm
    os_total = charge - ion_q
    if ve_total % n or os_total % n:
        raise ValueError(f"per-metal count is not a whole number ({ve_total}/{n} e, OS {os_total}/{n})")
    ve, ox = ve_total // n, os_total // n
    dn = group - ox
    if dn < 0:
        raise ValueError(f"d-electron count {dn} < 0 (oxidation state {ox} above group {group})")
    ve_ionic_total = n * dn + ion_e + 2 * mm
    if ve_ionic_total != ve_total:
        raise ValueError(f"covalent {ve_total} and ionic {ve_ionic_total} totals disagree")
    return {"metal": sym, "n": n, "group": group, "charge": charge, "mm": mm,
            "ve": ve, "os": ox, "dn": dn}


def fmt_os(ox):
    return "0" if ox == 0 else f"{'+' if ox > 0 else '−'}{abs(ox)}"
