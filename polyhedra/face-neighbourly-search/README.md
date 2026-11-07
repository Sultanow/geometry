# Face-neighbourly polyhedra: cluster search

Search for polyhedra in which every two faces share an edge (tetrahedron, Szilassi, Mizhaev, Röst–Vígh, …),
with eight and eleven faces. Every polyhedron whose faces lie in given planes is a choice of cells of the plane
arrangement; a SAT solver finds all choices that form a closed, embedded surface in which every two faces share an edge.

## Quality chain (nothing counts as a find before it has passed all of it)
1. SAT construction (`facenb/arrsat.py`): embedded by construction.
2. Floating-point check as in Röst–Vígh (`facenb/verify2.py`): faces simple; for every pair of faces, the intersection
   on the common line equals the shared edges.
3. Classification (`facenb/combtype.py`, `facenb/known.py`): exact isomorphism of the incidence structure against
   tetrahedron, Szilassi, Mizhaev, Röst–Vígh.
4. New types only: integer planes with the same type (`facenb/integerize.py`) and a proof in exact rational arithmetic
   (`facenb/exact.py`): general position, simple faces, closed, orientable, manifold, genus, embedded, face-neighbourly.

`tests/controls.py` must print `ALL CONTROLS PASSED` (reconstructs Mizhaev and Röst–Vígh from their planes, reproduces
Mizhaev's integer vertices and Röst–Vígh's exact volume, automorphism groups of order 4 and 8, rejects wrong geometries).

## Families (`facenb/families.py`)
C2_7 (control: must find Szilassi), S4_8 (Mizhaev's symmetry), D2_8 (Röst–Vígh's symmetry), C2_8, ASYM_8;
S4_11, D2_11 (4 + 4 + 2 + 1, slightly perturbed), C2_11 (1 + 5 pairs), ASYM_11.

## On rci-head
```bash
cd $HOME && unzip face-neighbourly-search.zip && cd face-neighbourly-search
bash slurm/setup.sh              # env + controls -> ALL CONTROLS PASSED
bash slurm/submit_A.sh           # eight faces: 384 tasks x 4 h
bash slurm/submit_B.sh           # eleven faces: 384 tasks x 12 h (queues behind A)
bash slurm/status.sh A           # queue + summary table (no exact step)
sbatch --export=ALL,RUN=A slurm/collect.sbatch   # final report: $FN_OUT/report_A/summary.md (+ new_type_<k>.json)
```
Paths in `slurm/config.sh`; logs in `$HOME/data/outputs/logs`.

## Guided search (run D)
Instead of yes/no, `facenb/guided.py` computes for a plane arrangement the largest number of pairs of faces that can share
an edge at once (MaxSAT; every face one simple polygon, encoded directly by reachability). `climb.py` moves the parameters of
a symmetric family uphill in this score (`facenb/gfamilies.py`); full scores go through the same checks as `search.py`.
On Mizhaev's family it needs about 16 times fewer evaluations per hit than blind sampling.
```bash
python tests/test_guided.py      # -> GUIDED CONTROLS PASSED
bash slurm/submit_D.sh           # 256 tasks x 12 h: 8 faces (calibration), 9, 10, 11 faces
bash slurm/progress.sh D         # best score per family out of all pairs, final scores of the climbs, full hits
bash slurm/status.sh D           # full hits classified as usual; final report with collect.sbatch (RUN=D)
```
