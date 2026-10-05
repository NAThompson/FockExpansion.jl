# Checks of the numerical ψ₄₀: closed-form Z⁰ and Z⁴ components, and the recurrence
# (Λ²-32)ψ₄₀ = 12ψ₄₁ + 2ψ₄₂ - 2Vψ₃₀ + 2Eψ₂₀ by finite differences.
using FockExpansion, Printf
const F=FockExpansion
E, a21 = -2.9, 0.4
p40_0(a,t)=(12*E^2-11*E+6*a21+1)/1152-(72*E*a21+E-30*a21-2)*sin(a)*cos(t)/720
p40_4(a,t)=sin(a)/36+7/288+2(4cos(a)^2-1)/(135π)
for (a,t) in ((0.7,1.1),(1.3,2.4))
    t0=@elapsed v0=psi40(a,t; Z=0.0, E, a21)
    @printf("(%.1f,%.1f) Z=0: %.15f closed %.15f diff %.1e (%.1f s)\n", a, t, v0, p40_0(a,t), v0-p40_0(a,t), t0)
    Zs=(0.5,1.0,1.5,2.0,2.5); vs=[psi40(a,t; Z, E, a21) for Z in Zs]
    A=[Z^k for Z in Zs, k in 0:4]; c=A\vs
    @printf("            Z⁴ coef: %.12f closed %.12f diff %.1e;  Z⁰ coef from fit %.12f\n", c[5], p40_4(a,t), c[5]-p40_4(a,t), c[1])
end
# recurrence by finite differences at Z=2
Z=2.0
ψ(a,t)=psi40(a,t; Z, E, a21, n=10, levels=10)
function rhs(a,t)
    ξ=F.xi_stable(a,t); V=1/ξ-Z*(1/sin(a/2)+1/cos(a/2))
    12psi41(a,t; Z, E, a21, n=10, levels=10)+2psi42(a,t; Z)-2V*psi30(a,t; Z, E, a21, rtol=1e-13)+2E*psi20(a,t; Z, E, a21)
end
for (a,t) in ((0.7,1.1),)
    h=2e-3; f0=ψ(a,t)
    faa=(ψ(a+h,t)-2*f0+ψ(a-h,t))/h^2; fa=(ψ(a+h,t)-ψ(a-h,t))/(2h)
    ftt=(ψ(a,t+h)-2*f0+ψ(a,t-h))/h^2; ft=(ψ(a,t+h)-ψ(a,t-h))/(2h)
    L=-4*(faa+2cot(a)*fa+(ftt+cot(t)*ft)/sin(a)^2)-32*f0
    @printf("FD recurrence at (%.1f,%.1f): (Λ²-32)ψ₄₀ = %.9f, rhs = %.9f, diff %.1e\n", a, t, L, rhs(a,t), L-rhs(a,t))
end
