"""Progress of guided runs: per family the best score reached (out of all n(n-1)/2 pairs), the histogram of final scores of
the climbs, evaluations, timeouts and full hits."""
import argparse, glob, json, os, re
ap = argparse.ArgumentParser(); ap.add_argument('--root', required=True); a = ap.parse_args()
fam = {}
for f in glob.glob(os.path.join(a.root, '*', 'stats_*.json')):
    s = json.load(open(f))
    if s.get('mode') != 'guided': continue
    t = fam.setdefault(s['family'], dict(tasks=0, cpu_h=0, evals=0, timeouts=0, climbs=0, best=-1, nmax=s.get('nmax'), hist={}, hits=0))
    t['tasks'] += 1; t['cpu_h'] += s['seconds'] / 3600; t['evals'] += s['samples']; t['timeouts'] += s['timeouts']; t['climbs'] += s['climbs']
    t['best'] = max(t['best'], s['best']); t['hits'] += s['solutions'] - s['float_rejected']; t['nmax'] = t['nmax'] or s.get('nmax')
    for k, v in s['final_hist'].items(): t['hist'][int(k)] = t['hist'].get(int(k), 0) + v
print('| family | tasks | CPU h | evaluations | timeouts | climbs | best / all pairs | final scores of climbs (top 5) | full hits |')
print('|---|---|---|---|---|---|---|---|---|')
for name, t in sorted(fam.items(), key=lambda kv: (int(re.findall(r'\d+', kv[0])[-1]), kv[0])):
    top = ', '.join(f'{k}: {v}' for k, v in sorted(t['hist'].items(), reverse=True)[:5])
    print(f"| {name} | {t['tasks']} | {t['cpu_h']:.1f} | {t['evals']} | {t['timeouts']} | {t['climbs']} | {t['best']} / {t['nmax']} | {top} | {t['hits']} |")
