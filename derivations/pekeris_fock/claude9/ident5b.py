import sys
sys.path.insert(0, '/Users/nathompson7/Documents/Codex/2026-09-05/how-much-compu/outputs/claude_work')
from pslq_grid import identify
from mpmath import mp, mpf, pi, catalan, log
mp.dps = 60
G, l2 = catalan, log(2)
B = {'G,log2': [mpf(1), 1/pi, 1/pi**2, G/pi, G/pi**2, l2/pi, l2/pi**2]}
for line in open('coef5.txt'):
    lab, val, err = [s.strip() for s in line.split('|')]
    if lab.startswith('1 ') and lab in ('1 r1^0 r2^0 r12^0', '1 r1^0 r2^0 r12^2', '1 r1^0 r2^0 r12^4', '1 r1^2 r2^0 r12^0'):
        hits = identify(mpf(val), 38, B, verbose=False, expansions=False)
        print(lab, [r for t, r in hits if t.startswith('x[x*1]')])
