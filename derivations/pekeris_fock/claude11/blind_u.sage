# Does the T5 reflection at u = -1 preserve the type-II image W1 at the circle point u = -sqrt(8/9) on x2 = 1/3?
from ore_algebra import *
R.<l> = QQ[]; A.<Dl> = OreAlgebra(R)
def op(fn): return sum(R(c) * Dl^m for m, c in enumerate(load(fn)))
LA = op('op_e2_Q.sobj')
aA = QQbar(sqrt(QQ(8)/9))
b0 = QQ(-97)/100
r = QQ(1)/200
def loop(L, base, a):
    s = sign(a - base); p0 = a - s*r
    return L.numerical_transition_matrix([base, p0, a - s*r*I, a + s*r, a + s*r*I, p0, base], 1e-40)
Tc = loop(LA, b0, -aA)
Tu = loop(LA, b0, -1)
I5 = identity_matrix(Tc.base_ring(), 5)
Nc = Tc - I5; Nu = Tu - I5
def mx(M): return max(abs(x) for x in M.list())
print('rank-ish of Tc-1: sing vals', [abs(x) for x in Nc.change_ring(CDF).SVD()[1].diagonal()])
print('Tu-1 sing vals', [abs(x) for x in Nu.change_ring(CDF).SVD()[1].diagonal()])
print('|(Tc-1)^2| =', mx(Nc*Nc))
print('|(Tu-1)(Tc-1)| =', mx(Nu*Nc), '  |(Tc-1)(Tu-1)| =', mx(Nc*Nu))
print('|Tu Tc - Tc Tu| =', mx(Tu*Tc - Tc*Tu))
# also the circle point +aA and u=+1 for comparison
b1 = QQ(97)/100
Tc2 = loop(LA, b1, aA); Tu2 = loop(LA, b1, 1)
print('+side: |(Tu-1)(Tc-1)| =', mx((Tu2-I5)*(Tc2-I5)), ' |[Tu,Tc]| =', mx(Tu2*Tc2 - Tc2*Tu2))
