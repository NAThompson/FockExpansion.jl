# I₁ = -(F(1) - F(c))/4 with F = ∫ T3(c) f(c) dc, split as
#   F(1) - F(c) = a1·log(1-ρ) + ac·log(1-ρc) + poly(c, ρ).
# Prints poly and ac as Julia code; a1 is entered in factored form by hand.
from sympy import *
c, r = symbols('c rho', positive=True)
L1, Lc = symbols('L1 Lc')
T3 = 4*c**3-3*c
for kind, f in (('inv', 1/(1-r*c)), ('log', log(1-r*c)), ('xlog', (1-r*c)*log(1-r*c))):
    F = expand(integrate(expand(T3*f), c))
    D = expand(F.subs(c, 1) - F)
    D = D.subs({log(r-1): L1, log(1-r): L1, log(c*r-1): Lc, log(1-c*r): Lc})
    D = expand(D)
    a1 = factor(D.coeff(L1)); ac = factor(D.coeff(Lc))
    poly = factor(expand(D - D.coeff(L1)*L1 - D.coeff(Lc)*Lc))
    assert not poly.has(log)
    print(f"# {kind}\n# a1 = {a1}")
    print(f"ac_{kind}(c, ρ) = {julia_code(ac).replace('rho','ρ')}")
    print(f"poly_{kind}(c, ρ) = {julia_code(poly).replace('rho','ρ')}")
