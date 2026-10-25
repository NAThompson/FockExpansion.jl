from ore_algebra import *
R.<l> = QQ[]; A.<Dl> = OreAlgebra(R)
cs = load('op_e2_Q.sobj')
L = sum(R(c) * Dl^m for m, c in enumerate(cs))
a = QQbar(sqrt(QQ(8)/9))
# loop around +sqrt(8/9) based at 1/5: go to a - 1/150, circle via 4 points, back
r = QQ(1)/150
pts = [QQ(1)/5, a - r, a - r*I, a + r, a + r*I, a - r, QQ(1)/5]
T = L.numerical_transition_matrix(pts, 1e-40)
print(T.base_ring(), T[0,0])
