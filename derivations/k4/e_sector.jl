# ∂_E ψ₄₀ from (H-E)∂_EΨ = Ψ:  (Λ²-32)∂_Eψ₄₀ = 12∂_Eψ₄₁ - 2V∂_Eψ₃₀ + 2E∂_Eψ₂₀ + 2ψ₂₀.
# The χ inside ψ₂₀ is handled by the resolvent identity: with (Λ²-12)χ = h_χ,
#   L₄⁺χ = (L₄⁺h_χ - χ + P₄χ)/20.   Checked here against differences of psi40.
using FockExpansion, Printf
const F=FockExpansion
Z, a21 = 2.0, 0.4
ξ(a,t)=F.xi_stable(a,t); σ(a)=cos(a/2)+sin(a/2)
function dEpsi40(α, θ, E; n=8, levels=8)
    function src(a, t)   # everything except 2Zχ
        s=sin(a); v=s*cos(t)
        V=1/ξ(a,t)-Z*(1/sin(a/2)+1/cos(a/2))
        φ41=Z*(π-2)*v/(30π)                              # ∂_E of Zψ⁽¹⁾
        φ30=Z*σ(a)*(2+sin(a)/2)/18-(6-ξ(a,t)^2)*ξ(a,t)/72  # r₁r₂ = sin α/2
        ψ20=(1-2E)/12+Z*(-σ(a)*ξ(a,t)/3)+Z^2*(1/3+s/2)+(a21-Z*(π+4)/(9π))*v
        12φ41-2V*φ30+2E*(-1/6)+2ψ20
    end
    hχ(a, t)=2σ(a)/(3ξ(a,t)*sin(a))-8(π-2)*sin(a)*cos(t)/(3π)
    base=F.solve_k4(src, α, θ; n, levels)
    Lh=F.solve_k4(hχ, α, θ; n, levels)
    base, Lh
end
pts=((0.7, 1.1), (1.3, 2.4), (0.4, 0.6), (1.1, 0.2))
rems=Float64[]; M=zeros(length(pts), 2)
for (i, (α, θ)) in enumerate(pts)
    E=-2.9; h=0.5
    fd=(psi40(α, θ; Z, E=E+h, a21)-psi40(α, θ; Z, E=E-h, a21))/(2h)   # ψ₄₀ is quadratic in E
    base, Lh=dEpsi40(α, θ, E)
    push!(rems, fd-(base+2Z*(Lh-F.chi(α, cos(θ)))/20))
    M[i, :]=[F.harmonic40(α, θ), F.harmonic42(α, θ)]
end
c=M[1:2, :]\rems[1:2]
@printf("P₄ part fitted on 2 points: c = (%.12f, %.12f); residual at the others: %s\n", c..., join([@sprintf("%.1e", r) for r in rems[3:end]-M[3:end, :]*c], ", "))
