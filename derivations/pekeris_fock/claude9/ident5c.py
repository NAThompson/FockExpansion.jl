import sys
sys.path.insert(0, '/Users/nathompson7/Documents/Codex/2026-09-05/how-much-compu/outputs/claude_work')
from pslq_grid import identify
from mpmath import mp, mpf, pi
mp.dps = 60
c = {l.split('|')[0].strip(): mpf(l.split('|')[1]) for l in open('coef5.txt')}
b = c['1 r1^0 r2^0 r12^4'] / mpf(1.5)
a = (c['1 r1^4 r2^0 r12^0'] - 2 * b) / 16
print('check r12^2 + 3b =', mp.nstr(c['1 r1^0 r2^0 r12^2'] + 3 * b, 5), ' r1^2 + 16a + 2b =', mp.nstr(c['1 r1^2 r2^0 r12^0'] + 16 * a + 2 * b, 5))
rest = c['1 r1^0 r2^0 r12^0'] - 3 * a - mpf(1.5) * b
print('rest', rest, 'a', a, 'b', b)
B = {'pi': [mpf(1), pi, pi**2, 1/pi, 1/pi**2]}
hits = identify(rest, 36, B, verbose=False, expansions=False)
print([r for t, r in hits if t.startswith('x[x*1]')])
