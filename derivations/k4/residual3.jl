# Exact residual of the truncated Fock series, from the recurrence:
# coefficient of R^(k-2) log^p R in (H-E)Ψ is ½[(Λ²-k(k+4))ψ_kp - h_kp] (missing ψ := 0).
# For Ψ⁽³⁾: (H-E)Ψ/R² = (Vψ30-Eψ20) + (Vψ31-Eψ21) log R - E R (ψ30 + ψ31 log R).
using FockExpansion, Printf
Z, E, a21 = 2.0, -2.9037243770341195983, 0.47674787900
α, θ = 2.2, 2.4
ξ=sqrt(1-sin(α)*cos(θ)); V=1/ξ-Z*(1/sin(α/2)+1/cos(α/2))
s0=V*psi30(α, θ; Z, E, a21, rtol=1e-13)-E*psi20(α, θ; Z, E, a21)
s1=V*psi31(α, θ; Z)-E*psi21(α, θ; Z)
@printf("s0 = %.6f, s1 = %.6f   (paper: 16.1067 + 1.26173 log R)\n", s0, s1)
p30=psi30(α, θ; Z, E, a21, rtol=1e-13); p31=psi31(α, θ; Z)
for (R, ref) in ((1e-6, 1.32479067459e-12), (1e-3, 7.3890323515e-06), (0.1, 0.128163060598))
    res=R^2*(s0+s1*log(R))-E*R^3*(p30+p31*log(R))
    @printf("R=%g  closed-form residual %.11e   paper (AD, 256-bit) %.11e\n", R, abs(res), ref)
end
