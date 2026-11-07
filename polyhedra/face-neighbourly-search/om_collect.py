"""Summary of the oriented-matroid tests: python om_collect.py --root DIR --total 203"""
import argparse, glob, json, os
ap = argparse.ArgumentParser(); ap.add_argument('--root', required=True); ap.add_argument('--total', type=int, default=203); a = ap.parse_args()
R = [json.load(open(f)) for f in glob.glob(os.path.join(a.root, 'result_*.json'))]
cnt = {k: sum(r['result'] == k for r in R) for k in ('SAT', 'UNSAT', 'undecided')}
lines = [f'# Oriented-matroid test: {len(R)}/{a.total} candidates done', '',
         f"- UNSAT (no realisation as an embedded polyhedron with planar faces): {cnt['UNSAT']}",
         f"- SAT (a chirotope satisfies all necessary conditions; realisation to be searched): {cnt['SAT']}  {sorted(r['candidate'] for r in R if r['result'] == 'SAT')}",
         f"- undecided within the time limit: {cnt['undecided']}  {sorted(r['candidate'] for r in R if r['result'] == 'undecided')}", '']
ts = sorted(r['seconds'] for r in R if r['result'] != 'undecided')
if ts: lines.append(f'- solving time of decided candidates: median {ts[len(ts) // 2]:.0f} s, max {ts[-1]:.0f} s')
open(os.path.join(a.root, 'summary.md'), 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
