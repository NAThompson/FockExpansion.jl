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

# The fourth-order coefficients below are computed in Float64. Refuse wider inputs rather
# than silently dropping their precision.
function check_float64(name, xs...)
    any(x->x isa AbstractFloat && precision(x)>precision(Float64), xs) &&
        throw(ArgumentError("$name is Float64-only; got $(join(unique(typeof.(xs)), ", "))"))
    nothing
end

# Unnormalized k = 4 hyperspherical harmonics with l = 0, 2 (l = 1 is odd under
# α → π-α and absent from singlet states).
harmonic40(α, θ) = 4cos(α)^2-1
harmonic42(α, θ) = sin(α)^2*(3cos(θ)^2-1)/2

# Coefficients of harmonic40 and harmonic42 in ψ₄₁, from the ψ₄₀ solvability
# condition ⟨Y₄ₗ, h₄₀⟩ = 0. The E and a₂₁ terms were identified exactly. The
# pure-Z terms involve the transcendental parts of ψ₃₀; they come from the 256-bit
# spectral solution of the hierarchy (derivations/k4/fock_spectral_ref.jl, 40 digits),
# rounded to Float64.
const PSI41_HARMONIC = (
    l0 = (-0.0013554801797822578866, 0.0030593204017684670484, 0.00070190494196250615738),
    l2 = (0.0051053947955489353034, -0.0053679564244078672069, -0.0019003487766684624263),
)
function psi41_harmonic(α, θ; Z, E, a21)
    π=typedpi(α)
    c0=evalpoly(Z, (0, PSI41_HARMONIC.l0...))
    c2=evalpoly(Z, (0, PSI41_HARMONIC.l2...))
    # E and a₂₁ enter only through Y₄₀ + 4Y₄₂ = 3(1-2sin²α sin²θ).
    c=Z*((3π-8)*E/(540π)-(5π-14)*a21/(90π))
    (c0+c)*harmonic40(α, θ)+(c2+4c)*harmonic42(α, θ)
end

"""Z² component of ψ₄₁ by quadrature against the S³ Green's function (about a second);
an independent check of `psi41_z2`."""
function psi41_z2_green(α::Real, θ::Real; n = 8, levels = 8)
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
a40 and a42 are not fixed by the Fock recurrence. Float64 only; one two-dimensional
tanh-sinh quadrature with ψ₃₀ at every node: `step = 0.0625` (default) gives about 1e-14
in about a second, `step = 0.125` about 1e-10 in 0.3 s. Inputs wider than Float64 throw
an ArgumentError."""
function psi40(α, θ; Z, E, a21, a40 = 0.0, a42 = 0.0, step = 0.0625)
    check_float64(:psi40, α, θ, Z, E, a21, a40, a42)
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
    solve_k4_symmetric(source, α, θ; h2, step)+a40*harmonic40(α, θ)+a42*harmonic42(α, θ)
end

"""ψ₄₁, the coefficient of R⁴log R, at an interior angle.
Includes its Y₄ₗ part, which is fixed by the ψ₄₀ equation and depends on E and a₂₁.
Float64 only (about 15 ms); inputs wider than Float64 throw an ArgumentError.
The Y₄ₗ constants are exact to Float64 rounding; the remaining error is that of the Z²
component ψ₄₁⁽²⁾ (about 1e-15), so the absolute error grows like Z²."""
function psi41(α, θ; Z, E, a21)
    check_float64(:psi41, α, θ, Z, E, a21)
    g=geometry(α, θ)
    (; α, θ, v, ξ)=g
    s=sin(α)
    z1=(π-2)/(2880π)*(3*(32*E-15)-8*(12*E-5)*ξ^2)
    z3=-(π-2)/(120π)*(4+5s)*v
    z2=iszero(Z) ? 0.0 : psi41_z2(α, θ)
    Z*z1+Z^2*z2+Z^3*z3+psi41_harmonic(α, θ; Z, E, a21)
end
