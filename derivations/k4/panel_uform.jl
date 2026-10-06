# The ψ₃₀ panel kernel Â/√P is a superposition of inverse linear forms on S³:
#   (2/π) atan(√P/(-2yt))/√P = (2/π) ∫₀¹ (-2yt) du / (4y²t² + P u²),
# and for b = cos α, y = -ξ/√2, 4y²t² + Pu² = C(t,u)(1 - q(t,u)·x) with |q| < 1.
using FockExpansion, QuadGK, Printf
const F=FockExpansion
for (b, y) in ((0.3, -0.6), (-0.4, -0.7))
    direct=first(F.panel_adaptive(b, y, 1.0, 0.0; rtol=1e-13))
    P(t)=1+b+(2-4y^2)*t^2+(1-b)*t^4
    h=sqrt((1+b)/(1-b)); upper=min(h, 1.0)
    # panel_adaptive integrates over q ∈ (0, min(h,1)); compare the N₀ moment via its integrand form
    N0(t)=log((1+t)/(1-t))*(2/π)*atan(sqrt(P(t))/(-2y*t))/sqrt(P(t))
    Nu(t)=log((1+t)/(1-t))*(2/π)*quadgk(u->(-2y*t)/(4y^2*t^2+P(t)*u^2), 0, 1; rtol=1e-14)[1]
    a=quadgk(N0, 0, 1; rtol=1e-12)[1]; c=quadgk(Nu, 0, 1; rtol=1e-12)[1]
    @printf("b=%.1f y=%.1f: ∫ L Â/√P = %.13f,  u-form %.13f\n", b, y, a, c)
end
