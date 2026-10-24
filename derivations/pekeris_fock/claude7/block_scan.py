import json, sys
sys.path.insert(0, '../claude6')
from ratrec import ratrec
xs = []; vals = None; meta = None
for k in range(8):
    d = json.load(open(f'bl_{k}.json')); xs += d['xs']; meta = d['meta']
    if vals is None: vals = {n: [] for n in d['vals']}
    for n in d['vals']: vals[n] += d['vals'][n]
order = sorted(range(len(xs)), key=lambda i: xs[i]); xs = [xs[i] for i in order]
print('points', len(xs))
for name in vals:
    V = [vals[name][i] for i in order]; keys, tg = meta[name]
    worst = (0, 0); unres = 0; nnz = 0
    for ti in range(len(tg)):
        for j in range(len(keys)):
            ys = [v[ti][j] for v in V]
            if not any(ys): continue
            nnz += 1
            r = ratrec(xs, ys)
            if r is None: unres += 1
            else: worst = (max(worst[0], r[0]), max(worst[1], r[1]))
    print(f'{name}: rank {len(keys)}, closure targets {len(tg)}, nonzero coefficients {nnz}, max line degree (num,den) {worst}, unresolved {unres}', flush=True)
