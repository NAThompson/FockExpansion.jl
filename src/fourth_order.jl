# Fourth-order logarithmic Fock coefficients ψ₄₂ (R⁴log²R) and ψ₄₁ (R⁴log R).
#
# Both solve (Λ²-32)ψ₄ₚ = h₄ₚ with
#   h₄₂ = 0,   h₄₁ = 24ψ₄₂ - 2Vψ₃₁ + 2Eψ₂₁.
# k = 4 is resonant: Λ² has eigenvalue 32 on the harmonics Y₄ₗ, so a regular
# solution exists only when h is orthogonal to them. That condition for h₄₁
# fixes ψ₄₂, and the same condition for h₄₀ = 12ψ₄₁ + 2ψ₄₂ - 2Vψ₃₀ + 2Eψ₂₀
# fixes the Y₄ₗ part of ψ₄₁. Only the Y₄ₗ part of ψ₄₀ is left free.
#
# By powers of Z, ψ₄₁ = Zψ⁽¹⁾ + Z²ψ⁽²⁾ + Z³ψ⁽³⁾ + harmonic part. ψ⁽¹⁾ and ψ⁽³⁾ are
# the closed forms of Liverts & Barnea, PRA 92, 042512 (2015), Table I; ψ⁽²⁾ is
# computed here by quadrature against the S³ Green's function.

"""ψ₄₂, the coefficient of R⁴log²R (Abbott–Maslen; Liverts–Barnea 2015, Eq. 206).
With r₁, r₂ the electron positions, 1-2sin²α sin²θ = 1-8|r₁×r₂|²/R⁴."""
function psi42(α, θ; Z)
    π=typedpi(α)
    g=geometry(α, θ)
    (π-2)*(5π-14)/(180π^2)*Z^2*(1-2sin(g.α)^2*sin(g.θ)^2)
end

# Unnormalized k = 4 hyperspherical harmonics with l = 0, 2 (l = 1 is odd under
# α → π-α and absent from singlet states).
harmonic40(α, θ) = 4cos(α)^2-1
harmonic42(α, θ) = sin(α)^2*(3cos(θ)^2-1)/2

# Coefficients of harmonic40 and harmonic42 in ψ₄₁, from the ψ₄₀ solvability
# condition ⟨Y₄ₗ, h₄₀⟩ = 0. The E and a₂₁ terms were identified exactly; the
# pure-Z terms involve the transcendental parts of ψ₃₀ and are numerical
# (derivations/k4/harmonic41.jl, accurate to about 1e-13).
const PSI41_HARMONIC = (
    l0 = (-0.001355480179782, 0.003059320401767, 0.000701904941965),
    l2 = (0.005105394795546, -0.005367956424392, -0.001900348776693),
)
function psi41_harmonic(α, θ; Z, E, a21)
    π=typedpi(α)
    c0=evalpoly(Z, (0, PSI41_HARMONIC.l0...))
    c2=evalpoly(Z, (0, PSI41_HARMONIC.l2...))
    # E and a₂₁ enter only through Y₄₀ + 4Y₄₂ = 3(1-2sin²α sin²θ).
    c=Z*((3π-8)*E/(540π)-(5π-14)*a21/(90π))
    (c0+c)*harmonic40(α, θ)+(c2+4c)*harmonic42(α, θ)
end

"""Z² component of ψ₄₁ without Y₄ₗ admixture, by quadrature against the S³ Green's function."""
function psi41_z2(α::Real, θ::Real; n = 8, levels = 8)
    solve_k4(psi41_z2_source, Float64(α), Float64(θ); n, levels)
end

# Source of ψ₄₁⁽²⁾ without its Y₄ part (which the Green's function removes anyway).
function psi41_z2_source(a, t)
    B=(π-2)/(3π)
    ξ=xi_stable(a, t)
    ς=cos(a/2)+sin(a/2)
    V1=1/sin(a/2)+1/cos(a/2)
    B*(V1*(5ξ^3/6-ξ)+ς*(ξ-1/ξ))
end

"""ψ₄₀, the R⁴ coefficient, at an interior angle: the solution with no Y₄₀, Y₄₂
component plus `a40*Y₄₀ + a42*Y₄₂` (Y₄₀ = 4cos²α-1, Y₄₂ = sin²α P₂(cos θ)).
a40 and a42 are not fixed by the Fock recurrence. Float64 only; one
two-dimensional quadrature with ψ₃₀ at every node (a few seconds)."""
function psi40(α, θ; Z, E, a21, a40 = 0.0, a42 = 0.0, n = 8, levels = 8)
    g=geometry(α, θ)
    α, θ=Float64(g.α), Float64(g.θ)
    # (Λ²-32)ψ₄₀ = 12ψ₄₁ + 2ψ₄₂ - 2Vψ₃₀ + 2Eψ₂₀. The harmonic parts of ψ₄₁ and ψ₄₂
    # are annihilated by the pure inverse; 12Z²ψ₄₁⁽²⁾ = 12Z²(Λ²-32)⁺ h₄₁⁽²⁾ is applied
    # through the iterated kernel.
    function source(a, t)
        ξ=xi_stable(a, t)
        s, c=sin(a), cos(t)
        V=1/ξ-Z*(1/sin(a/2)+1/cos(a/2))
        z1=(π-2)/(2880π)*(3*(32*E-15)-8*(12*E-5)*ξ^2)
        z3=-(π-2)/(120π)*(4+5s)*s*c
        12*(Z*z1+Z^3*z3)-2V*psi30(a, t; Z, E, a21, rtol = 1e-12)+2E*psi20(a, t; Z, E, a21)
    end
    h2=iszero(Z) ? nothing : (a, t)->12Z^2*psi41_z2_source(a, t)
    solve_k4(source, α, θ; n, levels, h2)+a40*harmonic40(α, θ)+a42*harmonic42(α, θ)
end

"""ψ₄₁, the coefficient of R⁴log R, at an interior angle.
Includes its Y₄ₗ part, which is fixed by the ψ₄₀ equation and depends on E and a₂₁.
Float64 only; each call takes a two-dimensional quadrature (about a second)."""
function psi41(α, θ; Z, E, a21, n = 8, levels = 8)
    g=geometry(α, θ)
    (; α, θ, v, ξ)=g
    s=sin(α)
    z1=(π-2)/(2880π)*(3*(32*E-15)-8*(12*E-5)*ξ^2)
    z3=-(π-2)/(120π)*(4+5s)*v
    z2=iszero(Z) ? 0.0 : psi41_z2(α, θ; n, levels)
    Z*z1+Z^2*z2+Z^3*z3+psi41_harmonic(α, θ; Z, E, a21)
end
