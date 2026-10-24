import json, sys
from guess_line import ratinterp, P
d = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'guess_line.json'))
xs, vals, targets = d['xs'], d['vals'], d['targets']
n = len(xs)
for ti, tg in enumerate(targets):
    nz = [j for j in range(45) if any(v[ti][j] for v in vals)]
    out = []
    for j in nz:
        ys = [v[ti][j] for v in vals]
        found = None
        for tot in range(0, n - 2):
            for dd in range(tot + 1):
                if ratinterp(xs, ys, tot - dd, dd): found = (tot - dd, dd); break
            if found: break
        out.append((j, found))
    print('target', tg, 'nonzero coefficients:', len(nz), ' degrees (num,den):', out, flush=True)
