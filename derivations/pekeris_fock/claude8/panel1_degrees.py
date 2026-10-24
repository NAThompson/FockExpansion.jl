import sys
sys.path.insert(0, '../claude6')
from ratrec import ratrec
from panel1_fp import sym_transform, P
xs = []; vals = []
lam = 0
while len(xs) < 120:
    lam += 1
    h, Y2 = (5 + lam) % P, (7 + 2 * lam) % P
    try:
        mats = [sym_transform(h, Y2, ax, fam) for ax in (0, 1) for fam in ('A', 'D')]
    except (ValueError, StopIteration):
        continue
    xs.append(lam); vals.append(mats)
nz = 0; worst = (0, 0); unres = 0
for m in range(4):
    for i in range(39):
        for j in range(39):
            ys = [v[m][i][j] for v in vals]
            if not any(ys): continue
            nz += 1
            r = ratrec(xs, ys)
            if r is None: unres += 1
            else: worst = (max(worst[0], r[0]), max(worst[1], r[1]))
print('points', len(xs), 'nonzero entries', nz, 'max line degree (num,den)', worst, 'unresolved', unres)
