import sys
sys.path.insert(0, '/Users/nathompson7/Documents/Codex/2026-09-05/how-much-compu/outputs/claude_work')
from pslq_grid import identify
from mpmath import mp, mpf, pi
mp.dps = 60
B = {'pi': [mpf(1), pi, pi**2, 1/pi, 1/pi**2]}
nz = idok = 0
out = open('coef5_ident.txt', 'w')
for line in open('coef5.txt'):
    lab, val, err = [s.strip() for s in line.split('|')]
    x = mpf(val)
    if abs(x) < mpf(10)**-35:
        continue
    nz += 1
    hits = identify(x, 36, B, verbose=False, expansions=False)
    h = [r for t, r in hits if t.startswith('x[x*1]')]
    if h:
        idok += 1
    out.write(f"{lab} | {mp.nstr(x, 20)} | {h[0] if h else 'NONE'}\n")
print('nonzero', nz, 'identified (plain x, scale 1)', idok)
