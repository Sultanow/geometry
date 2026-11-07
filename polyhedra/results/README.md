# Results

Outputs of the cluster runs behind the theorem that there is no face-neighbourly nine-faced polyhedron of genus 3.
The SAT formulas and DRAT proofs themselves (about 1.8 GB) are on Kaggle:
https://doi.org/10.34740/kaggle/dsv/21808309

| Path | Content |
|---|---|
| `enum_g3/candidates.json` | the 203 candidate maps. A map is a list of 9 faces; a face is the cycle of its vertices, each vertex the triple `[i, j, k]` of the faces meeting there. The list index is the candidate number used everywhere. |
| `enum_g3/cube_<k>.json`, `enum_g3/summary.md` | the complete enumeration (256 cubes, 108,096 labelled solutions, 203 classes) and its cross-check against the first enumeration |
| `om_cert_g3/cert_<k>.json` | decision for candidate k: UNSAT in the v2 encoding, DRAT proof verified by drat-trim, UNSAT in the alternative encoding, timings |
| `om_cert_g3/summary.md` | 203/203 UNSAT in both encodings, 203/203 proofs verified |
| `analysis_F/`, `analysis_G/` | near-misses with 34 of 36 adjacent pairs, exactly verified with integer planes (`analysis.jsonl`, `summary.md`) |
| `runs_A-G.tar.gz` | raw statistics, hits and best configurations of the geometric search runs A–G (blind sampling, guided search, growth from the known eight-faced polyhedra, second generation) |

Produced with the code in [`../face-neighbourly-search/`](../face-neighbourly-search/) on the Research Computing Infrastructure of the EDIH Hamburg.
