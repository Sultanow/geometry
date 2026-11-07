"""Summary of the certified run: python om_cert_collect.py --root DIR --total 203"""
import argparse, glob, json, os
ap = argparse.ArgumentParser(); ap.add_argument('--root', required=True); ap.add_argument('--total', type=int, default=203); a = ap.parse_args()
R = [json.load(open(f)) for f in glob.glob(os.path.join(a.root, 'cert_*.json'))]
unsat = [r for r in R if r['result'] == 'UNSAT']; ver = [r for r in unsat if r.get('drat_trim') == 'VERIFIED']
alt = [r for r in R if r.get('alt_result') == 'UNSAT']; both = [r for r in ver if r.get('alt_result') == 'UNSAT']
lines = [f'# Certified oriented-matroid test: {len(R)}/{a.total} candidates done', '',
         f'- UNSAT with the v2 encoding: {len(unsat)}',
         f'- of these, DRAT proof VERIFIED by drat-trim: {len(ver)}',
         f'- UNSAT with the alternative encoding (cross-check): {len(alt)}',
         f'- UNSAT in both encodings AND DRAT-verified: {len(both)}',
         f"- SAT (v2): {sorted(r['candidate'] for r in R if r['result'] == 'SAT')}",
         f"- disagreement between the encodings: {sorted(r['candidate'] for r in R if r['result'] != r.get('alt_result'))}",
         f"- proof not verified: {sorted(r['candidate'] for r in unsat if r.get('drat_trim') != 'VERIFIED')}"]
if len(both) == a.total: lines += ['', '**All candidates are excluded: every one is UNSAT in both encodings, with a verified DRAT proof.**']
open(os.path.join(a.root, 'summary.md'), 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
