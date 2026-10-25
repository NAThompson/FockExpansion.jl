# Evaluate the normalised telescoper coefficients c_m(lam0)/c_r(lam0) mod p at many lam0 (one job = a list of lam0 values).
import sys
NO_MAIN = True
load('gd.sage')
spec, r, extra, start, count, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), sys.argv[7]
with open(out, 'w') as f:
    for t in range(start, start + count):
        lam0 = QQ(t + 3) / 1  # integer sample points mod p
        res = solve_at(lam0, r, 3*r + extra, 8*r + 2*extra)
        if res is None or len(res[0]) != 1:
            f.write('%d FAIL\n' % (t + 3)); f.flush(); continue
        v = res[0][0]
        v = [x / v[r] for x in v] if v[r] != 0 else v
        f.write('%d %s\n' % (t + 3, ' '.join(str(x) for x in v))); f.flush()
