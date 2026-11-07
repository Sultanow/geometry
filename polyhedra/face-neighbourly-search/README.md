# face-neighbourly-search

Code behind the theorem that there is no face-neighbourly nine-faced polyhedron of genus 3 (every two faces share an
edge, three faces at every vertex). Data and certificates: https://doi.org/10.34740/kaggle/dsv/21808309

## Pipeline of the proof

| Step | Code | Cluster script | Output |
|---|---|---|---|
| 1. all combinatorial types (dual triangulations, 9 vertices, genus 3) | `facenb/enumerate2.py`, `enumerate_run.py`, `enumerate_collect.py` | `slurm/enum.sbatch`, `slurm/enum_collect.sbatch` | `../results/enum_g3/` (203 classes) |
| 1'. independent cross-check (K9 embeddings with three quadrilaterals) | `facenb/enumerate.py` | inside `enumerate_collect.py --crosscheck` | `../results/enum_g3/summary.md` |
| 2. oriented-matroid SAT test with DRAT proofs, second formulation as cross-check | `facenb/omsat.py`, `om_cert.py`, `om_cert_collect.py` | `slurm/setup_drat.sh`, `slurm/om_cert.sbatch` | `../results/om_cert_g3/` (203 UNSAT, 203 verified) |

`facenb/omsat.py` contains both formulations: `add_map` + `add_reflex` (v2) and `add_map_alt` + `add_reflex_alt`
(other ray lines). Conditions: Grassmann-Pluecker relations, simple faces, no edge piercing a face, at most one reflex
corner per vertex, reorientation symmetry breaking.

## Geometric searches (near-misses, background)

| Code | Purpose |
|---|---|
| `facenb/arrsat.py` | all face-neighbourly polyhedra whose faces lie in given planes (SAT on the cells of the arrangement) |
| `search.py`, `facenb/families.py` | blind sampling of symmetric plane families (runs A-C) |
| `facenb/guided.py`, `facenb/gfamilies.py`, `climb.py` | guided search: MaxSAT score, genus-aware and face-size-aware scores, growth from the known eight-faced polyhedra (runs D-G) |
| `analyze_best.py` | recompute, verify and analyse the best configurations (`../results/analysis_F`, `analysis_G`) |
| `facenb/realize.py`, `facenb/realize_fast.py`, `facenb/realize_om.py` | realisation of a fixed combinatorial type by penalty minimisation (CMA-ES), optionally guided by a chirotope. This is where the gap case was found that led to the reflex condition. |
| `facenb/verify2.py`, `facenb/exact.py`, `facenb/integerize.py` | independent float check (interval test) and exact verification in rational arithmetic with integer planes |
| `facenb/combtype.py`, `facenb/known.py` | combinatorial types; the known polyhedra (Mizhaev, Roest-Vigh, Szilassi) |

## Controls

Every result is only trusted after these pass:

```bash
python tests/controls.py          # ALL CONTROLS PASSED   (SAT construction, exact verification)
python tests/test_guided.py       # GUIDED CONTROLS PASSED
python tests/test_enumerate.py    # ENUMERATION CONTROLS PASSED
python tests/test_omsat.py        # OM CONTROLS PASSED
python tests/test_omsat_alt.py    # OM CONTROLS PASSED (alternative encoding)
```

Requirements: Python 3.12, numpy, shapely, networkx, python-sat; for the realisation modules also numba and cma.
The proof checker drat-trim is a git submodule in `tools/drat-trim` (`git clone --recurse-submodules`, then
`make -C tools/drat-trim`), or built by `bash slurm/setup_drat.sh`.

## Cluster

Paths are set in `slurm/config.sh` (Slurm, partition `std`, conda environment `facenb`). `slurm/setup.sh` creates the
environment and runs the controls.

## License

MIT (code). Data in `../results/` and on Kaggle: CC BY 4.0.
