# Iterated Green's function G2 = G∘G of Δ+8 on S³, zonal, orthogonal to degree-2 harmonics:
#   (Δ+8)G2 = G,  G(γ) = [(γ-π)cos3γ + sin(3γ)/6]/(4π² sinγ).  With G2 = w/sinγ: w'' + 9w = sinγ·G.
from sympy import *
g = symbols('gamma', real=True)
rhs = ((g-pi)*cos(3*g)+sin(3*g)/6)/(4*pi**2)
w = Function('w')
sol = dsolve(w(g).diff(g,2)+9*w(g)-rhs, w(g)).rhs
C1, C2 = symbols('C1 C2')
# regularity: w(0) = 0 and w(π) = 0
s1 = solve([sol.subs(g,0), sol.subs(g,pi)], [C1, C2], dict=True)
print("conditions:", s1)
wsol = simplify(sol.subs(s1[0]) if s1 else sol)
print("w =", wsol)
# remaining free multiple of sin3γ: fix by orthogonality to U2 = sin3γ/sinγ on S³ (weight sin²γ):
c = symbols('c')
wfull = wsol + c*sin(3*g)
free = [s for s in wfull.free_symbols if s not in (g,)]
print("free symbols:", free)
cond = integrate(expand(wfull*sin(3*g)), (g, 0, pi))
cs = solve(cond, [s for s in free][-1])
print("orthogonality:", cs)
