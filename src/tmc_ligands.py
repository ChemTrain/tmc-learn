"""Ligand library for TMC Learn.

Each ligand has a class in Green's covalent bond classification (CBC):
L donates two electrons, X one, so "L2X" (Cp) gives 2*2 + 1 = 5 electrons
when counted as a neutral ligand. Everything else follows from the class:

    covalent (neutral) count   2L + X electrons, no charge
    ionic count                2L + 2X electrons, charge -X
    oxidation state            each X adds +1 to the metal

so Cl (X) gives 1 e covalent and 2 e as Cl⁻, Cp (L2X) gives 5 e or 6 e as Cp⁻,
and an oxo (X2) gives 2 e or 4 e as O²⁻.

The one exception is linear NO, written `ionic=(2, +1)`: covalent it is a
3-electron LX ligand, and in the ionic count it is NO⁺, a 2-electron donor
that lowers the metal's oxidation state by one. That is the convention of
most textbooks; complexes with NO are flagged so that their oxidation state
is never asked (see tmc_sets.py).

Fields
    name     what the ligand is called on the page
    show     short formula or abbreviation, as in a formula (PPh₃, Cp*)
    cls      CBC class: L, X, LX, L2, X2, L2X, L3, LX2, X3, L2X2
    smiles   a drawing fragment; the atoms that bind the metal are mapped
             :1, :2, ... Atoms mapped :51 and up are drawn as text labels
             (see `labels`), which keeps PPh₃ or Mes from filling the page.
    points   how the mapped atoms bind: one tuple per bond to the metal.
             (1,) a σ bond from atom 1; (1, 2) an η²-bond to the midpoint of
             atoms 1 and 2; (1, 2, 3, 4, 5) an η⁵-ring, drawn to its centre.
             Default: every mapped atom is its own σ bond.
    order    bond order drawn to the metal at each point (1 dative/single,
             2 for M=CR₂, M=O, M=NR, 3 for M≡N, M≡CR). Default all 1.
    labels   {map number: (label when the ligand points right, label when it
             points left)}; RDKit markup <sub></sub> is allowed.
    bridging True for μ-ligands, which bind two metals. Their class is the
             class of the whole ligand (μ-Cl is LX: 2 e to one metal, 1 to
             the other in the neutral count).
    ionic    (electrons, charge) when the ionic count does not follow from
             the class (only linear NO).
    ion      the ligand's name as an ion in the ionic count, when that reads
             better than the formula (Cl⁻, Cp⁻, O²⁻).
    xy, mpos hand-set 2D coordinates (one per atom, bond length about 1.5)
             and where the metal sits in them, for ligands RDKit cannot lay
             out bound (an s-cis diene, a tub-shaped cod, an η³-allyl).
             Chelates bound through single atoms need neither: build.py lays
             them out with a stand-in metal bonded to their donors.
    depictor "rdkit" to lay the ligand out with RDKit's own depictor instead
             of CoordGen, when that draws it better (salen).
"""
import re

PH3 = ("PPh<sub>3</sub>", "Ph<sub>3</sub>P")


def lab(right, left=None):
    return (right, left or right)


LIGANDS = {
    # ---------------- L: two-electron neutral donors ----------------
    "CO":     dict(name="carbonyl", show="CO", cls="L", smiles="[C:1]#O", labels={1: lab("C")}),
    "CNtBu":  dict(name="tert-butyl isocyanide", show="CNt-Bu", cls="L", smiles="[C:1]#N[*:51]",
                   labels={1: lab("C"), 51: lab("tBu")}),
    "PPh3":   dict(name="triphenylphosphine", show="PPh₃", cls="L", smiles="[P:1]", labels={1: PH3}),
    "PMe3":   dict(name="trimethylphosphine", show="PMe₃", cls="L", smiles="[P:1]",
                   labels={1: ("PMe<sub>3</sub>", "Me<sub>3</sub>P")}),
    "PEt3":   dict(name="triethylphosphine", show="PEt₃", cls="L", smiles="[P:1]",
                   labels={1: ("PEt<sub>3</sub>", "Et<sub>3</sub>P")}),
    "PCy3":   dict(name="tricyclohexylphosphine", show="PCy₃", cls="L", smiles="[P:1]",
                   labels={1: ("PCy<sub>3</sub>", "Cy<sub>3</sub>P")}),
    "PiPr3":  dict(name="triisopropylphosphine", show="P(i-Pr)₃", cls="L", smiles="[P:1]",
                   labels={1: ("P(iPr)<sub>3</sub>", "(iPr)<sub>3</sub>P")}),
    "PtBu3":  dict(name="tri-tert-butylphosphine", show="P(t-Bu)₃", cls="L", smiles="[P:1]",
                   labels={1: ("P(tBu)<sub>3</sub>", "(tBu)<sub>3</sub>P")}),
    "POPh3":  dict(name="triphenyl phosphite", show="P(OPh)₃", cls="L", smiles="[P:1]",
                   labels={1: ("P(OPh)<sub>3</sub>", "(PhO)<sub>3</sub>P")}),
    "NH3":    dict(name="ammine", show="NH₃", cls="L", smiles="[NH3:1]"),
    "H2O":    dict(name="aqua", show="H₂O", cls="L", smiles="[OH2:1]"),
    "py":     dict(name="pyridine", show="py", cls="L", smiles="[n:1]1ccccc1"),
    "Clpy":   dict(name="3-chloropyridine", show="3-Cl-py", cls="L", smiles="[n:1]1cc(Cl)ccc1"),
    "MeCN":   dict(name="acetonitrile", show="MeCN", cls="L", smiles="[N:1]#CC"),
    "THF":    dict(name="tetrahydrofuran", show="THF", cls="L", smiles="[O:1]1CCCC1"),
    "N2":     dict(name="dinitrogen", show="N₂", cls="L", smiles="[N:1]#N"),
    "H2":     dict(name="dihydrogen (η²)", show="η²-H₂", cls="L", smiles="[H:1][H:2]", points=[(1, 2)],
                   xy=[(-0.45, 0), (0.45, 0)], mpos=(0, -1.5)),
    "C2H4":   dict(name="ethylene (η²)", show="η²-C₂H₄", cls="L", smiles="[CH2:1]=[CH2:2]", points=[(1, 2)]),
    "PhCCPh": dict(name="diphenylacetylene (η², 2 e)", show="η²-PhC≡CPh", cls="L",
                   smiles="c1ccccc1[C:1]#[C:2]c1ccccc1", points=[(1, 2)]),
    "SIMes":  dict(name="SIMes (saturated N-heterocyclic carbene)", show="SIMes", cls="L",
                   smiles="[C:1]1N([*:51])CCN1[*:52]", labels={51: lab("Mes"), 52: lab("Mes")}),
    "IMes":   dict(name="IMes (N-heterocyclic carbene)", show="IMes", cls="L",
                   smiles="[C:1]1N([*:51])C=CN1[*:52]", labels={51: lab("Mes"), 52: lab("Mes")}),
    "IPr":    dict(name="IPr (N-heterocyclic carbene)", show="IPr", cls="L",
                   smiles="[C:1]1N([*:51])C=CN1[*:52]", labels={51: lab("Dipp"), 52: lab("Dipp")}),
    "CHPh":   dict(name="benzylidene (Fischer-type carbene in Grubbs catalysts)", show="=CHPh", cls="L",
                   smiles="[CH:1]c1ccccc1", order=[2]),
    "COMePh": dict(name="methoxy(phenyl)carbene (Fischer carbene)", show="=C(OMe)Ph", cls="L",
                   smiles="[C:1](OC)c1ccccc1", order=[2]),
    "COMeMe": dict(name="methoxy(methyl)carbene (Fischer carbene)", show="=C(OMe)Me", cls="L",
                   smiles="[C:1](OC)C", order=[2]),
    "OPPh3":  dict(name="triphenylphosphine oxide", show="OPPh₃", cls="L", smiles="[O:1]=[P:51]",
                   labels={51: PH3}),
    # chelating and π-bonded neutral ligands
    "dppe":   dict(name="dppe (1,2-bis(diphenylphosphino)ethane)", show="dppe", cls="L2",
                   smiles="[P:1]CC[P:2]", labels={1: ("PPh<sub>2</sub>", "Ph<sub>2</sub>P"), 2: ("PPh<sub>2</sub>", "Ph<sub>2</sub>P")}),
    "bpy":    dict(name="2,2′-bipyridine", show="bpy", cls="L2", smiles="c1cc[n:1]c(c1)-c1cccc[n:2]1"),
    "phen":   dict(name="1,10-phenanthroline", show="phen", cls="L2", smiles="c1c[n:1]c2c(c1)ccc1ccc[n:2]c12"),
    "en":     dict(name="ethylenediamine", show="en", cls="L2", smiles="[NH2:1]CC[NH2:2]"),
    "dach":   dict(name="1,2-diaminocyclohexane", show="dach", cls="L2", smiles="[NH2:1]C1CCCCC1[NH2:2]"),
    "terpy":  dict(name="2,2′:6′,2″-terpyridine", show="terpy", cls="L3",
                   smiles="c1cc[n:1]c(c1)-c1cccc(-c2cccc[n:3]2)[n:2]1"),
    "cod":    dict(name="1,5-cyclooctadiene (η²,η²)", show="cod", cls="L2",
                   smiles="[CH:1]1=[CH:2]CC[CH:3]=[CH:4]CC1", points=[(1, 2), (3, 4)],
                   # a tub seen from the side: both C=C face the metal
                   xy=[(-0.55, 0), (-1.35, 0.55), (-0.75, 1.9), (0.75, 1.9), (1.35, 0.55), (0.55, 0), (0.75, 1.2), (-0.75, 1.2)],
                   mpos=(0, -1.0)),
    "nbd":    dict(name="norbornadiene (η²,η²)", show="nbd", cls="L2",
                   smiles="[CH:1]1=[CH:2]C2[CH:3]=[CH:4]C1C2", points=[(1, 2), (3, 4)]),
    "C4H6":   dict(name="1,3-butadiene (η⁴)", show="η⁴-C₄H₆", cls="L2",
                   smiles="[CH2:1]=[CH:2][CH:3]=[CH2:4]", points=[(1, 2, 3, 4)],
                   xy=[(-1.5, 0), (-0.75, 1.3), (0.75, 1.3), (1.5, 0)], mpos=(0, -1.4)),
    "C4H4":   dict(name="cyclobutadiene (η⁴)", show="η⁴-C₄H₄", cls="L2",
                   smiles="[CH:1]1=[CH:2][CH:3]=[CH:4]1", points=[(1, 2, 3, 4)]),
    "C6H6":   dict(name="benzene (η⁶)", show="η⁶-C₆H₆", cls="L3",
                   smiles="[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1", points=[(1, 2, 3, 4, 5, 6)]),
    "cymene": dict(name="p-cymene (η⁶)", show="η⁶-p-cymene", cls="L3",
                   smiles="[c:1]1(C)[cH:2][cH:3][c:4](C(C)C)[cH:5][cH:6]1", points=[(1, 2, 3, 4, 5, 6)]),
    "C6Me6":  dict(name="hexamethylbenzene (η⁶)", show="η⁶-C₆Me₆", cls="L3",
                   smiles="[c:1]1(C)[c:2](C)[c:3](C)[c:4](C)[c:5](C)[c:6]1C", points=[(1, 2, 3, 4, 5, 6)]),
    "HGcarb": dict(name="2-isopropoxybenzylidene (chelating, Hoveyda–Grubbs)", show="=CH(C₆H₄-o-Oi-Pr)", cls="L2",
                   smiles="[CH:1]c1ccccc1[O:2]C(C)C", order=[2, 1]),

    # ---------------- X: one-electron (anionic) ligands ----------------
    "H":      dict(name="hydride", show="H", cls="X", smiles="[H:1]", ion="H⁻"),
    "F":      dict(name="fluoride", show="F", cls="X", smiles="[F:1]", ion="F⁻"),
    "Cl":     dict(name="chloride", show="Cl", cls="X", smiles="[Cl:1]", ion="Cl⁻"),
    "Br":     dict(name="bromide", show="Br", cls="X", smiles="[Br:1]", ion="Br⁻"),
    "I":      dict(name="iodide", show="I", cls="X", smiles="[I:1]", ion="I⁻"),
    "Me":     dict(name="methyl", show="Me", cls="X", smiles="[CH3:1]", ion="CH₃⁻",
                   labels={1: ("CH<sub>3</sub>", "H<sub>3</sub>C")}),
    "Et":     dict(name="ethyl", show="Et", cls="X", smiles="[CH2:1]C", ion="C₂H₅⁻"),
    "CH2tBu": dict(name="neopentyl", show="CH₂t-Bu", cls="X", smiles="[CH2:1]C(C)(C)C", ion="neopentyl⁻"),
    "Ph":     dict(name="phenyl", show="Ph", cls="X", smiles="[c:1]1ccccc1", ion="Ph⁻"),
    "Cy":     dict(name="cyclohexyl", show="Cy", cls="X", smiles="[CH:1]1CCCCC1", ion="Cy⁻"),
    "CF3":    dict(name="trifluoromethyl", show="CF₃", cls="X", smiles="[C:1](F)(F)F", ion="CF₃⁻",
                   labels={1: lab("C")}),
    "COMe":   dict(name="acetyl", show="C(O)Me", cls="X", smiles="[C:1](C)=O", ion="acetyl⁻"),
    "CN":     dict(name="cyanide", show="CN", cls="X", smiles="[C:1]#N", ion="CN⁻", labels={1: lab("C")}),
    "OH":     dict(name="hydroxide", show="OH", cls="X", smiles="[OH:1]", ion="OH⁻"),
    "OMe":    dict(name="methoxide", show="OMe", cls="X", smiles="[O:1]C", ion="MeO⁻"),
    "OiPr":   dict(name="isopropoxide", show="Oi-Pr", cls="X", smiles="[O:1]C(C)C", ion="i-PrO⁻"),
    "OtBu":   dict(name="tert-butoxide", show="Ot-Bu", cls="X", smiles="[O:1]C(C)(C)C", ion="t-BuO⁻"),
    "ORF6":   dict(name="hexafluoro-tert-butoxide", show="OCMe(CF₃)₂", cls="X",
                   smiles="[O:1][*:51]",
                   labels={51: ("CMe(CF<sub>3</sub>)<sub>2</sub>", "(F<sub>3</sub>C)<sub>2</sub>MeC")}, ion="(CF₃)₂MeCO⁻"),
    "OAc":    dict(name="acetate (κ¹)", show="OAc", cls="X", smiles="[O:1]C(C)=O", ion="AcO⁻"),
    "NMe2":   dict(name="dimethylamide", show="NMe₂", cls="X", smiles="[N:1](C)C", ion="Me₂N⁻"),
    "NObent": dict(name="nitrosyl, bent", show="NO (bent)", cls="X", smiles="[N:1]=O", ion="NO⁻",
                   xy=[(0, 0), (1.05, 1.07)], mpos=(0, -1.5)),

    # ---------------- X2, X3: multiply bonded and dianionic chelates ----------------
    "O":      dict(name="oxo", show="O", cls="X2", smiles="[O:1]", order=[2], ion="O²⁻"),
    "CH2":    dict(name="methylidene (Schrock alkylidene)", show="=CH₂", cls="X2", smiles="[CH2:1]", order=[2], ion="CH₂²⁻",
                   labels={1: ("CH<sub>2</sub>", "H<sub>2</sub>C")}),
    "CHtBu":  dict(name="neopentylidene (Schrock alkylidene)", show="=CHt-Bu", cls="X2",
                   smiles="[CH:1]C(C)(C)C", order=[2], ion="t-BuCH²⁻"),
    "CHCMe2Ph": dict(name="neophylidene (Schrock alkylidene)", show="=CHCMe₂Ph", cls="X2",
                     smiles="[CH:1]C(C)(C)c1ccccc1", order=[2], ion="PhMe₂CCH²⁻"),
    "ox":     dict(name="oxalate (κ²)", show="ox", cls="X2", smiles="[O:1]C(=O)C(=O)[O:2]", ion="C₂O₄²⁻"),
    "cbdca":  dict(name="cyclobutane-1,1-dicarboxylate (κ²)", show="cbdca", cls="X2",
                   smiles="[O:1]C(=O)C1(CCC1)C(=O)[O:2]", ion="cbdca²⁻"),
    "N":      dict(name="nitrido", show="N", cls="X3", smiles="[N:1]", order=[3], ion="N³⁻"),
    "CtBu":   dict(name="neopentylidyne (alkylidyne)", show="≡Ct-Bu", cls="X3", smiles="[C:1]C(C)(C)C",
                   order=[3], ion="t-BuC³⁻"),

    # ---------------- LX, L2X, LX2, L2X2 ----------------
    "allyl":  dict(name="allyl (η³)", show="η³-C₃H₅", cls="LX", smiles="[CH2:1]=[CH:2][CH2:3]",
                   points=[(1, 2, 3)], ion="C₃H₅⁻", xy=[(-1.25, 0), (0, 0.75), (1.25, 0)], mpos=(0, -1.55)),
    "acac":   dict(name="acetylacetonate (κ²)", show="acac", cls="LX", smiles="[O:1]=C(C)C=C(C)[O:2]", ion="acac⁻"),
    "OAc2":   dict(name="acetate (κ²)", show="κ²-OAc", cls="LX", smiles="[O:1]=C(C)[O:2]", ion="AcO⁻"),
    "NO":     dict(name="nitrosyl, linear", show="NO", cls="LX", smiles="[N:1]#O", ionic=(2, +1), ion="NO⁺"),
    "TsDPEN": dict(name="TsDPEN (N-tosyl-1,2-diphenylethylenediamide)", show="TsDPEN", cls="LX",
                   smiles="[NH2:1]C([*:52])C([*:53])[N:2][*:51]", labels={51: lab("Ts"), 52: lab("Ph"), 53: lab("Ph")},
                   ion="TsDPEN⁻"),
    "Cp":     dict(name="cyclopentadienyl (η⁵)", show="Cp", cls="L2X",
                   smiles="[cH:1]1[cH:2][cH:3][cH:4][cH:5]1", points=[(1, 2, 3, 4, 5)], ion="Cp⁻"),
    "Cpstar": dict(name="pentamethylcyclopentadienyl (η⁵)", show="Cp*", cls="L2X",
                   smiles="[c:1]1(C)[c:2](C)[c:3](C)[c:4](C)[c:5]1C", points=[(1, 2, 3, 4, 5)], ion="Cp*⁻"),
    # no drawing yet: a κ³ tripod needs its own layout
    "Tp":     dict(name="hydrotris(pyrazolyl)borate (κ³)", show="Tp", cls="L2X", smiles=None, ion="Tp⁻"),
    "NAr":    dict(name="arylimido, linear (2,6-diisopropylphenyl)", show="=NAr", cls="LX2",
                   smiles="[N:1]c1c(C(C)C)cccc1C(C)C", order=[2], ion="ArN²⁻"),
    "NtBu":   dict(name="tert-butylimido, linear", show="=Nt-Bu", cls="LX2", smiles="[N:1]C(C)(C)C",
                   order=[2], ion="t-BuN²⁻"),
    "salen":  dict(name="salen (κ⁴-N,N,O,O)", show="salen", cls="L2X2",
                   smiles="[O:1]c1ccccc1C=[N:2]CC[N:3]=Cc1ccccc1[O:4]", ion="salen²⁻", depictor="rdkit"),
    "jsalen": dict(name="Jacobsen's salen (κ⁴-N,N,O,O)", show="salen*", cls="L2X2",
                   smiles="[O:1]c1c(C(C)(C)C)cc(C(C)(C)C)cc1C=[N:2]C1CCCCC1[N:3]=Cc1cc(C(C)(C)C)cc(C(C)(C)C)c1[O:4]",
                   ion="salen*²⁻"),
    "por":    dict(name="porphyrinato (κ⁴)", show="por", cls="L2X2",
                   smiles="c1cc2cc3ccc(cc4ccc(cc5ccc(cc1[n:1]2)[n:2]5)[n:3]4)[n:4]3", ion="por²⁻"),

    # ---------------- bridging ----------------
    "mu-Cl":  dict(name="chloride, bridging (μ-Cl)", show="μ-Cl", cls="LX", smiles="[Cl:1]", bridging=True, ion="Cl⁻"),
    "mu-Br":  dict(name="bromide, bridging (μ-Br)", show="μ-Br", cls="LX", smiles="[Br:1]", bridging=True, ion="Br⁻"),
    "mu-I":   dict(name="iodide, bridging (μ-I)", show="μ-I", cls="LX", smiles="[I:1]", bridging=True, ion="I⁻"),
    "mu-CO":  dict(name="carbonyl, bridging (μ-CO)", show="μ-CO", cls="L", smiles="[C:1]=O", bridging=True,
                   labels={1: lab("C")}),
    "mu-PhCCPh": dict(name="diphenylacetylene, bridging (μ-η²:η², 4 e)", show="μ-PhC≡CPh", cls="L2",
                      smiles="c1ccccc1[C:1]#[C:2]c1ccccc1", points=[(1, 2)], bridging=True),
    "mu-H":   dict(name="hydride, bridging (μ-H)", show="μ-H", cls="X", smiles="[H:1]", bridging=True, ion="H⁻"),
}

def parse_class(cls):
    """'L2X' -> (2, 1)."""
    if not re.fullmatch(r"(L\d?)?(X\d?)?", cls) or not cls:
        raise ValueError(f"bad CBC class {cls!r}")
    l = re.search(r"L(\d?)", cls)
    x = re.search(r"X(\d?)", cls)
    nl = int(l.group(1) or 1) if l else 0
    nx = int(x.group(1) or 1) if x else 0
    return nl, nx


def contribution(key):
    """(covalent electrons, ionic electrons, ionic charge, X count) for one ligand."""
    lig = LIGANDS[key]
    nl, nx = parse_class(lig["cls"])
    cov = 2 * nl + nx
    ion_e, ion_q = lig.get("ionic") or (2 * nl + 2 * nx, -nx)
    return cov, ion_e, ion_q, nx
