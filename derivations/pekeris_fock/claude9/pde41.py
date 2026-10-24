# Check that the closed form satisfies the psi41^(2) equation: Delta_6 (R^4 psi) = -R^2 [B(V1(5xi^3/6 - xi) + vs(xi - 1/xi)) + 24 psi42^(2)],
# with symbolic derivatives (Clausen: d/dz Cl2(z) = -log(2 sin(z/2))) evaluated at 100 digits.
import pickle, sympy as S
from mpmath import mp
mp.dps = 100
r1, r2, r12, pi = [S.Symbol(n, positive=True) for n in ('r1', 'r2', 'r12', 'pi')]
a_h, b_h = S.symbols('a_h b_h')
class Cl2(S.Function):
    def fdiff(self, argindex=1):
        return -S.log(2 * S.sin(self.args[0] / 2))
T = pickle.load(open('closed41.pkl', 'rb'))
T = {k: v.subs({S.Symbol('r1', positive=True): r1}) for k, v in T.items()}
R2 = r1**2 + r2**2; R = S.sqrt(R2); P = S.sqrt(2 * R2 - r12**2); d, s = r1 - r2, r1 + r2; L5 = r12 * P + s * d
al = 2 * S.atan(r2 / r1); be = S.asin((R2 - r12**2) / R2)
DL = Cl2(al - be) + Cl2(pi + be - al) - Cl2(pi - al - be) - Cl2(al + be)
lg = be * S.log((2 * r1 * r2)**2 / L5**2) + (be + al) * S.log((s**2 - r12**2) / (2 * r1 * r2)) + (be - al) * S.log((r12**2 - d**2) / (2 * r1 * r2))
A = {'1': 1, 'alpha': al, 'L': S.log((s + r12) / R), 'U': lg - DL + pi * S.log((d + r12)**2 / L5), 'V': P * (be + pi / 2)}
F = sum(T[k] * A[k] for k in T)            # = R^4 psi41^(2) (homogeneous degree 4)
D = lambda f, x: S.diff(f, x)
lap = (D(D(F, r1), r1) + 2 / r1 * D(F, r1) + D(D(F, r2), r2) + 2 / r2 * D(F, r2) + 2 * (D(D(F, r12), r12) + 2 / r12 * D(F, r12))
       + (r1**2 - r2**2 + r12**2) / (r1 * r12) * D(D(F, r1), r12) + (r2**2 - r1**2 + r12**2) / (r2 * r12) * D(D(F, r2), r12))
B = (pi - 2) / (3 * pi)
xi = r12 / R; V1 = R / r1 + R / r2 ; vs = (r1 + r2) / R   # 2sin(a/2)=|x-e1| etc: sin(a/2)=r2/R, cos(a/2)=r1/R
src = B * (V1 * (5 * xi**3 / 6 - xi) + vs * (xi - 1 / xi))
sa2st2 = 1 - ((r1**2 - r2**2) / R2)**2 - ((R2 - r12**2) / R2)**2     # sin^2 a sin^2 t = 1 - cos^2 a - (sin a cos t)^2
psi42 = (pi - 2) * (5 * pi - 14) / (180 * pi**2) * (1 - 2 * sa2st2)
expr = lap + R2 * (src + 24 * psi42)
f = S.lambdify((r1, r2, r12, pi, a_h, b_h), expr, modules=[{'Cl2': lambda z: mp.clsin(2, z)}, 'mpmath'])
import random
random.seed(3)
for _ in range(4):
    x1, x2 = mp.mpf(random.uniform(0.3, 1.2)), mp.mpf(random.uniform(0.3, 1.2))
    x12 = abs(x1 - x2) + (x1 + x2 - abs(x1 - x2)) * mp.mpf(random.uniform(0.1, 0.9))
    print(mp.nstr(f(x1, x2, x12, mp.pi, mp.mpf('0.123'), mp.mpf('-0.456')), 5), ' (scale', mp.nstr(abs(F.subs({}) if False else 1), 2), ')')
