<img src="logo.svg" alt="ChemTrain logo" height="110" align="right">

# TMC Learn

TMC Learn is a Quizlet-style trainer for counting electrons in transition-metal complexes. Everything is in one HTML file.

**Play it here:** https://chemtrain.github.io/tmc-learn/tmc_learn.html

You can also download `tmc_learn.html` and open it in any browser.

## What's inside

- 199 complexes in six sets: textbook essentials, metal carbonyls, Cp/arene/allyl complexes, 16- and 14-electron complexes, hydrides and metal–carbon multiple bonds, and Werner complexes; neutral, cationic and anionic, from 8 to 22 electrons
- A Ligands set: 47 ligands with their L/X class and the electrons they give, each drawn on a generic metal to show how it binds
- A 2D structure for nearly every complex, with the metal highlighted
- Two counting methods, covalent (neutral ligands, Green's L/X classification) and ionic; every answer shows the tally by both
- Learn (multiple choice, then typed answers) asks for the electron count, the oxidation state or dⁿ; flashcards; space for your own cards or a Quizlet export
- Geometry: from the number of ligands, the electron count and dⁿ, guess the shape (linear to octahedral) of 133 complexes; on Hard the look-alike with the same number of ligands is always among the choices
- A fold-out d-block periodic table with group numbers, the metal of the current question highlighted
- Sort: deal 6, 9 or 12 complexes and put each in its box: under 18, 18 or over 18 electrons (Easy), or 14, 15, 16, 17, 18, 19 (Hard)
- Count it: work through the count of a complex line by line, one ligand at a time, then the metal, the charge and the total (or the oxidation state and dⁿ in the ionic count)
- Name it: 47 complexes with a common name (Wilkinson's, Vaska's, Grubbs, Schrock, Zeise's salt, …); name the structure or pick the structure for the name, then a bonus question on its electron count
- Easy and Hard levels: answer choices are at least 2 apart on Easy and next to each other on Hard; in Name it, Hard takes the wrong answers from the same family

## Counting conventions

- Linear NO counts as a 3-electron LX ligand (NO⁺ in the ionic count). Its complexes are never asked for an oxidation state.
- Grubbs alkylidenes count as neutral 2-electron carbenes (Ru(II)); Schrock alkylidenes, oxo and imido count as X₂ or LX₂. Linear imido gives 4 electrons; π-donation from oxo, alkoxide and amide is not counted.
- Cyclobutadiene counts as a neutral 4-electron L₂ ligand.
- A bridging CO gives 2 electrons shared by two metals; each metal–metal bond gives one electron to each metal. Complexes with several metals are counted per metal.

## Sources

Every count is worked out by both methods from the complex's ligands, and the build stops if it differs from the textbook value typed in the set. The counts and conventions follow R. H. Crabtree, *The Organometallic Chemistry of the Transition Metals*; J. F. Hartwig, *Organotransition Metal Chemistry*; and M. L. H. Green and G. Parkin, [*J. Chem. Educ.* 2014, 91, 807](https://doi.org/10.1021/ed400504f).

Structures are drawn with [RDKit](https://www.rdkit.org).

## How to cite

If you use TMC Learn in a course, a worksheet or a paper, please cite it:

> Sala, F. *TMC Learn* (2026). ChemTrain. https://github.com/ChemTrain/tmc-learn

GitHub's **Cite this repository** button (from `CITATION.cff`) gives the same reference in APA and BibTeX.

## License

- Code (the HTML, CSS and JavaScript): [MIT](LICENSE).
- Content (the card sets, ligand library, questions, notes and structure drawings): [CC BY 4.0](LICENSE-CONTENT.md). You may reuse and adapt it, as long as you credit TMC Learn.
- Logo (`logo.svg`): the locomotive is adapted from [Twemoji](https://github.com/jdecked/twemoji), © Twitter, Inc and other contributors, licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); the benzene-ring smoke was added for ChemTrain.

TMC Learn is a sibling of [pKa Learn](https://github.com/ChemTrain/pka-learn) and shares its page design.
