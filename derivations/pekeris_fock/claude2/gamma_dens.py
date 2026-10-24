import sympy as sp
import order2_atoms as oa
G = oa.gamma_table()
dens = set()
for a in oa.ATOMS:
    for k in range(3):
        for t, c in G[a][k].items():
            dens.add(sp.factor(sp.fraction(sp.cancel(c))[1]))
D = sp.Integer(1)
for d in dens:
    D = sp.lcm(D, d)
print('Gamma denominators:', dens)
print('lcm:', sp.factor(D))
