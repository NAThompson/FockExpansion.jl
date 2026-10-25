import os
p = 2147483549; os.environ['PRIME'] = str(p)
NO_MAIN = True
import sys
for spec, fn, lam in (("e2:x1=l,x2=1/3", 'op_e2_Q.sobj', QQ(2)/7), ("e2:x1=1/5,x2=l", 'op_eC_Q.sobj', QQ(3)/7)):
    sys.argv = ['gd.sage', spec, '5']
    load('gd.sage')
    res = solve_at(lam, 5, 15, 40)
    v = res[0][0]; k = GF(p)
    cs = load(fn)
    w = [k(c(lam)) for c in cs]
    print(spec, 'fresh prime match:', [x / v[5] for x in v] == [x / w[5] for x in w])
