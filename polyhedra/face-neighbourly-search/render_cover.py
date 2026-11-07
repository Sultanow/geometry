"""Cover image for the Kaggle dataset: the near-miss polyhedron (nine faces, 34 of 36 pairs adjacent).
Usage:  python render_cover.py analysis_G/analysis.jsonl cover.png
Needs numpy and matplotlib. Uses the exactly verified integer planes of the first configuration that has them."""
import sys, json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
src, out = sys.argv[1], (sys.argv[2] if len(sys.argv) > 2 else 'cover.png')
recs = [json.loads(l) for l in open(src) if l.strip()]
rec = next((r for r in recs if r.get('exact')), next(r for r in recs if r.get('planes')))
P = np.array(rec['exact']['planes'] if rec.get('exact') else rec['planes'], float)
faces = [[tuple(t) for t in f] for f in rec['faces']]
X = {t: np.linalg.solve(P[list(t), :3], P[list(t), 3]) for f in faces for t in f}
pts = np.array(list(X.values())); c = pts.mean(0); r = np.percentile(np.linalg.norm(pts - c, axis=1), 90)
cols = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#1f78b4', '#b2df8a']
fig = plt.figure(figsize=(12.8, 6.4), dpi=100, facecolor='#0e1117')
for k, (el, az) in enumerate([(24, 30), (24, 150)]):
    ax = fig.add_axes([0.02 + 0.49 * k, 0.04, 0.47, 0.86], projection='3d'); ax.set_facecolor('#0e1117')
    ax.add_collection3d(Poly3DCollection([[X[t] for t in f] for f in faces], facecolors=[cols[p % 9] for p in range(len(faces))],
                                         edgecolors='#f0f0f0', linewidths=0.5, alpha=0.88))
    ax.set_xlim(c[0] - r, c[0] + r); ax.set_ylim(c[1] - r, c[1] + r); ax.set_zlim(c[2] - r, c[2] + r)
    ax.view_init(el, az); ax.set_axis_off(); ax.set_box_aspect((1, 1, 1), zoom=1.6)
fig.text(0.5, 0.93, 'Nine planar faces, 34 of 36 pairs adjacent \u2014 and no nine-hedron of genus 3 with all 36',
         ha='center', color='#f0f0f0', fontsize=17)
plt.savefig(out, facecolor=fig.get_facecolor()); print('written', out)
