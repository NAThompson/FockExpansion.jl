# ψ₃₀ depends on a₂₁ only through -a₂₁ψ₃₁/(ZB), B = (π-2)/(3π). Feeding that into
# h₄₀ and using h₄₁ = 24ψ₄₂ - 2Vψ₃₁ + 2Eψ₂₁ gives  ∂ψ₄₀/∂a₂₁ = -(ψ₄₁ - harmonic)/(ZB).
using FockExpansion, Printf
const F=FockExpansion
B=(π-2)/(3π); Z, E = 2.0, -2.9
for (α, θ) in ((0.7, 1.1), (1.3, 2.4))
    d=psi40(α, θ; Z, E, a21=0.5)-psi40(α, θ; Z, E, a21=0.0)
    pure=psi41(α, θ; Z, E, a21=0.0)-F.psi41_harmonic(α, θ; Z, E, a21=0.0)
    @printf("(%.1f,%.1f): Δψ₄₀ = % .12f   -0.5·ψ₄₁ᵖᵘʳᵉ/(ZB) = % .12f   diff %.1e\n", α, θ, d, -0.5pure/(Z*B), d+0.5pure/(Z*B))
end
