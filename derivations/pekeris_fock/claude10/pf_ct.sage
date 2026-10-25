# Picard-Fuchs (telescoper) of the holomorphic period  ∫∫ dq dX / Y,  Y^2 = X^3 + 2B X^2 + A C X  (Codex x19 surface),
# on a target line, by two creative-telescoping steps (ore_algebra ct): first X, then q.
import sys, time
from ore_algebra import *
which = sys.argv[1]            # 'e2' or 'e1'
R.<X,q,l> = QQ[]
D = 1 + q^4; P1, P2 = q^4 - 1, 2*q^2
if which == 'e2':
    x1, x2 = l, QQ(1)/3
    ax, aP = x2, P2
else:
    x1, x2 = QQ(1)/5, l
    ax, aP = x1, P1
Px = P1*x1 + P2*x2
A = 1 - ax^2; B = aP - ax*Px; C = D^2 - Px^2
F = X^3 + 2*B*X^2 + A*C*X
Alg.<DX,Dq,Dl> = OreAlgebra(R)
I = Alg.ideal([2*F*DX + F.derivative(X), 2*F*Dq + F.derivative(q), 2*F*Dl + F.derivative(l)])
t0 = time.time()
T1 = I.ct(DX, certificates=False, infolevel=0)
print('step X done', time.time()-t0, [ (t.order() if hasattr(t,'order') else '?') for t in T1]); sys.stdout.flush()
R2.<q,l> = QQ[]
Alg2.<Dq,Dl> = OreAlgebra(R2)
J = Alg2.ideal([Alg2(str(t)) for t in T1])
T2 = J.ct(Dq, certificates=False, infolevel=0)
print('step q done', time.time()-t0)
for t in T2:
    print('TELESCOPER order', t.order(), 'degree', t.degree())
    save(t, 'pf_%s.sobj' % which)
    print(t)
