import json, sys
from ratrec import ratrec
d = json.load(open(sys.argv[1]))
xs, vals, targets = d['xs'], d['vals'], d['targets']
print('points', len(xs))
for ti, tg in enumerate(targets):
    nz = [j for j in range(45) if any(v[ti][j] for v in vals)]
    degs = []
    for j in nz:
        r = ratrec(xs, [v[ti][j] for v in vals])
        degs.append((j, None if r is None else (r[0], r[1])))
    print('target', tg, 'nonzero', len(nz), 'max deg (num,den):',
          (max((x[1][0] for x in degs if x[1]), default=None), max((x[1][1] for x in degs if x[1]), default=None)),
          'unresolved', sum(1 for x in degs if x[1] is None), flush=True)
