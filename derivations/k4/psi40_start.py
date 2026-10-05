# ψ₄₀ by powers of Z: (Λ²-32)ψ₄₀⁽ʲ⁾ = h₄₀⁽ʲ⁾,  h₄₀ = 12ψ₄₁ + 2ψ₄₂ - 2Vψ₃₀ + 2Eψ₂₀,
# V = V₀ - Z V₁, V₀ = 1/ξ, V₁ = 2ς/sin α. The Z⁰ and Z⁴ sources depend on one
# variable only and have elementary solutions; checked here symbolically.
from sympy import *
a, v, E, A21 = symbols('a v E a21', real=True)
s = sqrt(1-a**2); xi = sqrt(1-v); vs = sqrt(1+s)
def Lam2(f):
    return -4*((1-a**2)*diff(f,a,2)+(1-v**2)*diff(f,v,2)-2*a*v*diff(f,a,v)-3*a*diff(f,a)-3*v*diff(f,v))
V0, V1 = 1/xi, 2*vs/s
# Z⁰: ψ₃₀⁽⁰⁾ and ψ₂₀⁽⁰⁾ (a₂₁ = total coefficient of sin α cos θ in ψ₂₀)
p30_0 = xi*((1-2*E)/24+A21/2)+xi**3*((E-2)/72-5*A21/12)
p20_0 = (1-2*E)/12+A21*v
h0 = expand(simplify(-2*V0*p30_0+2*E*p20_0))
c0, c1 = symbols('c0 c1')
sol = solve(Poly(expand(Lam2(c0+c1*v)-32*(c0+c1*v)-h0), v).all_coeffs(), [c0, c1])
p40_0 = factor(sol[c0])+factor(sol[c1])*v
print("psi40^(0) =", p40_0)
print("  check:", simplify(Lam2(p40_0)-32*p40_0-h0))
# Z⁴: ψ₃₀⁽³⁾ = -ς(2+5 sin α)/36 (Liverts–Barnea Table I), the only Z⁴ source
p30_3 = -vs*(2+5*s)/36
h4 = 2*V1*p30_3
p40_4 = s/36+Rational(7,288)
print("psi40^(4) particular = sin(a)/36 + 7/288;  check:", simplify(Lam2(p40_4)-32*p40_4-h4))
# Y₄₀ projection of the particular solution (measure sin²α dα; θ integrates out)
al = symbols('alpha')
proj = integrate((sin(al)/36)*(4*cos(al)**2-1)*sin(al)**2, (al, 0, pi))/integrate((4*cos(al)**2-1)**2*sin(al)**2, (al, 0, pi))
print("pure psi40^(4) = sin(a)/36 + 7/288 -", simplify(proj), "* (4cos²a-1)")
