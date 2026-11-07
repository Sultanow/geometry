# Face-neighbourly polyhedra

A polyhedron is *face-neighbourly* if every two of its faces share an edge. Known examples: the tetrahedron, the Szilassi polyhedron (7 faces, genus 1), and two eight-faced polyhedra of genus 3 (Mizhaev 2026; Röst and Vígh 2026).

## Result

> **Theorem (computer-assisted).** There is no closed embedded polyhedron of genus 3 with nine planar faces in which every two faces share at least one edge and every vertex lies on exactly three faces.

- **All candidates:** exactly 203 combinatorial types, from two independent enumerations that agree.
- **None realisable:** for every type, an oriented-matroid SAT formula (the signs of the 4×4 determinants of the nine face planes and the plane at infinity) is unsatisfiable. This holds in two independent formulations, and every proof is a DRAT certificate checked by drat-trim.
- **Near-miss:** an exactly verified nine-faced polyhedron of genus 4 in which 34 of the 36 pairs of faces share an edge.

Not covered: genus 4 and higher, and vertices with four or more faces.

**Status:** not yet peer-reviewed.

## Data and certificates

All formulas, DRAT proofs, the 203 candidate maps and the near-miss configurations (about 1.8 GB):
**https://doi.org/10.34740/kaggle/dsv/21808309** (Kaggle, DOI 10.34740/kaggle/dsv/21808309)

The dataset README explains every file and how to check a certificate with drat-trim.

## Code

[`face-neighbourly-search/`](face-neighbourly-search/) contains the complete code:

| Module | Purpose |
|---|---|
| `facenb/enumerate.py`, `facenb/enumerate2.py` | the two enumerations of the combinatorial types |
| `facenb/omsat.py` | oriented-matroid SAT encodings (two formulations) |
| `facenb/arrsat.py` | all face-neighbourly polyhedra on a given plane arrangement (SAT) |
| `facenb/exact.py` | exact verification in rational arithmetic |
| `facenb/guided.py`, `climb.py` | the guided geometric searches that produced the near-misses |
| `tests/` | controls that must pass before any result is trusted |
| `slurm/` | cluster scripts for the full runs |

```bash
cd face-neighbourly-search
pip install numpy shapely networkx python-sat
python tests/test_enumerate.py      # ENUMERATION CONTROLS PASSED
python tests/test_omsat.py          # OM CONTROLS PASSED
python tests/test_omsat_alt.py      # OM CONTROLS PASSED (alternative encoding)
```

## Acknowledgments

The computations were carried out on the Research Computing Infrastructure of the EDIH Hamburg. Development of the code and the analysis was assisted by Claude (Anthropic).

## Citation

Eldar Sultanow (2026). *No Face-Neighbourly Nine-Hedron of Genus Three* [Data set]. Kaggle. https://doi.org/10.34740/kaggle/dsv/21808309

## License

Code: MIT. Data on Kaggle: CC BY 4.0.
