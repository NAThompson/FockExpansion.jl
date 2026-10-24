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
# the closed forms of Liverts & Barnea, PRA 92, 042512 (2015), Table I; ψ⁽²⁾ is the
# closed form `psi41_z2` (derivations/pekeris_fock/THEORY.md, Phase 9), checked
# against the S³ Green's function and a one-dimensional Feynman integral.

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

# Y₄ₗ coefficients that make psi41_z2 orthogonal to Y₄₀, Y₄₂ (256-bit quadrature, 48 digits;
# derivations/pekeris_fock/codex/phase9/x19_report.md). Numerical, not identified.
const PSI41_Z2_PURE=(big"0.00082106684957556469415831715509798526043967205959",
                     big"0.01992734793080898549051527717103411701281763036468")

"""Z² component of ψ₄₁ without Y₄ₗ admixture, in closed form. With R = 1, r₁ = cos(α/2),
r₂ = sin(α/2), r₁₂ = ξ, v = sin α cos θ, β = asin v and P = √(1+v), it is a combination of
1, α, log(r₁+r₂+r₁₂), P(β+π/2) and a Clausen–log atom U with quartic polynomial coefficients
over ℚ(π), plus the numerical Y₄ₗ constants `PSI41_Z2_PURE`
(derivations/pekeris_fock/THEORY.md, Phase 9). Generic in the number type (Float64,
BigFloat, ForwardDiff duals); about 1e-16 absolute in Float64, also near the coalescences.
Liverts & Barnea (PRA 92, 042512 (2015); arXiv:1601.06964) give ψ₄₁⁽²ᵈ⁾ only as an infinite
Legendre series with explicit terms; this finite form is from derivations/pekeris_fock (phase 9)."""
function psi41_z2(α, θ)
    π=typedpi(α)
    (; α, θ, v, ξ, σ, opv)=geometry(fold_half(α, π), θ)   # ψ₄₁⁽²⁾ is exchange symmetric
    s, c=sincos(α)
    μ=cos(θ)
    P=sqrt(opv)
    G=ξ*P                                     # cos β, without cancellation near v = ±1
    β=abs(primal(μ))==1 ? μ*α : atan(v, G)
    d=sqrt(2*one(π))*sin(π/4-α/2)             # r₁ - r₂
    U=clausen_logs(α, μ, β, s, c, G)+π*log((d+ξ)^2/(G+c))
    k=(π-2)/π
    p=(5π^2-3π-14)/(360π^2)+k*(ξ*σ*(36+40s-60ξ^2)-20s*v-28s*c)/360+
      oftype(π, PSI41_Z2_PURE[1])*harmonic40(α, θ)+oftype(π, PSI41_Z2_PURE[2])*harmonic42(α, θ)
    p+7k/(45π)*s*c*α+k/18*(3-2s^2-4ξ^2+2ξ^4)*log(σ+ξ)-7k/(90π)*c*v*U-7k/(45π)*ξ*v*P*(β+π/2)
end

"""ψ₄₀, the R⁴ coefficient, for 0 < α < π and 0 ≤ θ ≤ π: the solution with no Y₄₀, Y₄₂
component plus `a40*Y₄₀ + a42*Y₄₂` (Y₄₀ = 4cos²α-1, Y₄₂ = sin²α P₂(cos θ)).
a40 and a42 are not fixed by the Fock recurrence. Float64 only; one two-dimensional
tanh-sinh quadrature with ψ₃₀ at every node: `step = 0.0625` (default) gives about 1e-14
in about a second, `step = 0.125` about 1e-10 in 0.3 s. Inputs wider than Float64 throw
an ArgumentError."""
function psi40(α, θ; Z, E, a21, a40 = 0.0, a42 = 0.0, step = 0.0625)
    check_float64(:psi40, α, θ, Z, E, a21, a40, a42)
    # θ = 0, π are coordinate poles, not boundaries; the source nodes stay interior.
    0<=θ<=π || throw(DomainError(θ, "psi40 requires 0 ≤ θ ≤ π."))
    α, θ=Float64(geometry(α, π/2).α), Float64(θ)
    # (Λ²-32)ψ₄₀ = 12ψ₄₁ + 2ψ₄₂ - 2Vψ₃₀ + 2Eψ₂₀. The harmonic parts of ψ₄₁ and ψ₄₂
    # are annihilated by the pure inverse; 12Z²ψ₄₁⁽²⁾ = 12Z²(Λ²-32)⁺ h₄₁⁽²⁾ is applied
    # through the iterated kernel.
    function source(a, t)
        ξ=xi_stable(a, t)
        s, c=sin(a), cos(t)
        V=1/ξ-Z*(1/sin(a/2)+1/cos(a/2))
        z1=(π-2)/(2880π)*(3*(32*E-15)-8*(12*E-5)*ξ^2)
        z3=-(π-2)/(120π)*(4+5s)*s*c
        12*(Z*z1+Z^3*z3)-2V*psi30(a, t; Z, E, a21)+2E*psi20(a, t; Z, E, a21)
    end
    h2=iszero(Z) ? nothing : (a, t)->12Z^2*psi41_z2_source(a, t)
    solve_k4_symmetric(source, α, θ; h2, step)+a40*harmonic40(α, θ)+a42*harmonic42(α, θ)
end

"""ψ₄₁, the coefficient of R⁴log R, at an interior angle.
Includes its Y₄ₗ part, which is fixed by the ψ₄₀ equation and depends on E and a₂₁.
Closed form (under a microsecond); supports ForwardDiff duals. Wider inputs throw an
ArgumentError, because the Y₄ₗ constants are stored only to Float64 rounding."""
function psi41(α, θ; Z, E, a21)
    check_float64(:psi41, α, θ, Z, E, a21)
    g=geometry(α, θ)
    (; α, θ, v, ξ)=g
    s=sin(α)
    z1=(π-2)/(2880π)*(3*(32*E-15)-8*(12*E-5)*ξ^2)
    z3=-(π-2)/(120π)*(4+5s)*v
    z2=psi41_z2(α, θ)
    Z*z1+Z^2*z2+Z^3*z3+psi41_harmonic(α, θ; Z, E, a21)
end
