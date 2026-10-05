# Fourth-order Fock coefficients ψ₄₂, ψ₄₁ (and a start on ψ₄₀)

Recurrence (Liverts & Barnea 2015, Eq. 103): `(Λ² - k(k+4))ψ_kp = h_kp`,
`h_kp = 2(k+2)(p+1)ψ_k,p+1 + (p+1)(p+2)ψ_k,p+2 - 2Vψ_k-1,p + 2Eψ_k-2,p`.

## Resonance and what it fixes

For k = 4, Λ² - 32 annihilates the harmonics Y₄ₗ, so a regular ψ₄ₚ exists only
if h₄ₚ ⊥ Y₄ₗ (l = 0, 2 for singlets; l = 1 is odd under α → π-α).

* `h₄₁ ⊥ Y₄ₗ` fixes ψ₄₂ (a pure Y₄ combination). `solv42.jl` confirms the
  published ψ₄₂ to 14 digits.
* `h₄₀ = 12ψ₄₁ + 2ψ₄₂ - 2Vψ₃₀ + 2Eψ₂₀ ⊥ Y₄ₗ` fixes the Y₄ₗ part of ψ₄₁.
  Liverts & Barnea present only the "pure" ψ₄₁ and call its harmonic part
  undetermined; it is not. It depends on Z, E, a₂₁ (`harmonic41.jl`):

      c₀ = -0.001355480179782 Z + 0.003059320401767 Z² + 0.000701904941965 Z³
      c₂ =  0.005105394795546 Z - 0.005367956424392 Z² - 0.001900348776693 Z³
      plus Z[(3π-8)E/(540π) - (5π-14)a₂₁/(90π)] (Y₄₀ + 4Y₄₂),
      Y₄₀ = 4cos²α-1, Y₄₂ = sin²α P₂(cosθ), Y₄₀ + 4Y₄₂ = 3(1 - 2sin²α sin²θ).

  The E and a₂₁ coefficients are exact (PSLQ, then confirmed to 13 digits);
  the pure-Z ones involve ψ₃₀'s transcendental parts.
* Only the Y₄ₗ part of ψ₄₀ is free (set by the global wavefunction).

## ψ₄₂ and ψ₄₁

    ψ₄₂ = (π-2)(5π-14)/(180π²) Z² (1 - 2sin²α sin²θ) = … Z² (1 - 8|r₁×r₂|²/R⁴)
    ψ₄₁ = Zψ⁽¹⁾ + Z²ψ⁽²⁾ + Z³ψ⁽³⁾ + c₀Y₄₀ + c₂Y₄₂
    ψ⁽¹⁾ = (π-2)/(2880π) [3(32E-15) - 8(12E-5)ξ²]          (LB Table I; check_table.py)
    ψ⁽³⁾ = -(π-2)/(120π) (4 + 5 sinα) sinα cosθ             (LB Table I; check_table.py)
    (Λ²-32)ψ⁽²⁾ = B[V₁(5ξ³/6 - ξ) + ς(ξ - 1/ξ)] + 24ψ₄₂⁽²⁾,   B = (π-2)/(3π)

ψ⁽²⁾ is the hard part. Liverts & Barnea give it as subcomponents 2b, 2c (closed)
and 2d (Legendre series tabulated to l = 10); their 2018 paper reports that the
2b and l = 0, 2 terms of 2d disagree with the Green's-function calculation.
Here ψ⁽²⁾ is computed directly with the S³ Green's function (`green4.jl`,
`src/green_k4.jl`): with x = (cos α, sin α cos θ, sin α sin θ e^{iφ}) ∈ S³,
Λ² = -4Δ_{S³}, and the Green's function of Δ+8 orthogonal to the degree-2
harmonics is G(γ) = [(γ-π)cos3γ + sin(3γ)/6]/(4π² sin γ). It reproduces ψ⁽¹⁾
and ψ⁽³⁾ to 1e-17, and finite differences confirm the ψ⁽²⁾ equation to 1e-7.

On S³ the Coulomb factors are chord distances to three points:
ξ = |x-e₂|/√2, 2sin(α/2) = |x-e₁|, 2cos(α/2) = |x+e₁|, so V₁ = 2/|x-e₁| + 2/|x+e₁|.

### ψ⁽²⁾ as a one-dimensional integral (`feynman41.jl`, `src/feynman_k4.jl`)

The source is B Σ± Q±(x)/(d± ξ) with Q± = -1/3 - (7/3)w + (5/3)w² ± uw
(u = cos α, w = sin α cos θ). Three steps make it elementary:

1. Feynman: 1/(dξ) = (1/(√2π)) ∫₀¹ dt/(√(t(1-t)) ℓ), ℓ = 1 - p·x, p± = (±t, 1-t, 0, 0):
   a point source moving along the arc from e₂ to ±e₁.
2. Polynomials as p-derivatives: x_i/ℓ = -∂ᵢ log ℓ and x_i x_j/ℓ = ∂ᵢ∂ⱼ(ℓ log ℓ).
3. 1/ℓ, log ℓ and ℓ log ℓ depend on x only through c = p̂·x = cos γ. With U = W/sin γ,
   W'' + 9W = -f sin γ/4 + κ sin 3γ, solved by variation of parameters. Only
   polynomials in cos γ, log ℓ and A(γ) = (2/√(1-ρ²)) atan((1+ρ)tan(γ/2)/√(1-ρ²)) occur.

So ψ⁽²⁾ = (B/(√2π)) ∫₀¹ dt/√(t(1-t)) Σ± [-U_{1/ℓ}/3 + (7/3)∂₂U_{log ℓ} + (5/3)∂₂²U_{ℓlog ℓ}
± ∂₁∂₂U_{ℓlog ℓ}] - 0.0048977716325816 Y₄₀ + 0.0145608238240513 Y₄₂, the
last two terms removing the harmonic part (fitted against the Green's function to
8e-16; no closed form found by PSLQ). Coefficients that vanish as ρ → 1 are kept in
factored form (`zonal_coeffs.py`), which makes the integrand stable in Float64 down to
t = 1e-32 (checked against 768-bit evaluation); it grows like log t at both ends.
A 608-node graded rule gives about 2e-15 in about 15 ms, versus 1.4 s for the
two-dimensional Green's-function quadrature. Each zonal solution is checked against
quadrature of its defining integrals (1e-16).

## Residual of the truncated series (`residual3.jl`, `residual_study.jl`)

The coefficient of R^(k-2) log^p R in (H-E)Ψ is ½[(Λ²-k(k+4))ψ_kp - h_kp],
with missing coefficients set to zero, so the residual of a truncation is
finite and exact. It reproduces the 256-bit automatic-differentiation values in
`paper/fock_30_coefficient.tex` to 12 digits.

    Ψ⁽³⁾:        R²[s₀ + s₁ log R] + O(R³ log R),  s₀ = Vψ₃₀-Eψ₂₀, s₁ = Vψ₃₁-Eψ₂₁
    +ψ₄₂:        R²[(s₀-ψ₄₂) + (s₁-12ψ₄₂) log R] + …   (removes the Y₄ part of s₁)
    +ψ₄₂+ψ₄₁:    R²[s₀ - 6ψ₄₁ - ψ₄₂] + O(R³ log²R)    (no log at order R²)
    +ψ₄₀:        O(R³ log²R)                           (any Y₄ part of ψ₄₀)

`residual_study.out` (Z = 2, ground-state E and a₂₁, a₄₀ = a₄₂ = 0), L² norm of
(H-E)Ψ/R² over the angles:

| R | Ψ⁽³⁾ | +ψ₄₂ | +ψ₄₂+ψ₄₁ | all of k = 4 |
|---|---:|---:|---:|---:|
| 1e-4 | 23.42 | 23.42 | 21.48 | 0.0017 |
| 1e-2 | 21.73 | 21.72 | 21.41 | 0.148 |
| 0.1 | 20.78 | 20.78 | 20.81 | 1.40 |
| 0.5 | 18.06 | 18.06 | 18.11 | 6.23 |

The log terms alone barely help: the R² residual is dominated by its non-log part
s₀ ≈ 16-18, which only ψ₄₀ removes.

## ψ₄₀ (`psi40_start.py`, `check_psi40.jl`, `src/fourth_order.jl`)

Numerically, ψ₄₀ = (Λ²-32)⁺[12ψ₄₁ + 2ψ₄₂ - 2Vψ₃₀ + 2Eψ₂₀] + a₄₀Y₄₀ + a₄₂Y₄₂. The
harmonic parts of ψ₄₁, ψ₄₂ drop out, and the Z²ψ₄₁ term uses the iterated kernel
G₂ = G∘G = [(γ²-2πγ)/(48π²) + 1/72 + 1/(864π²)] sin3γ/sinγ (`g2.py`), so one
two-dimensional quadrature suffices. Checks: the Z⁰ and Z⁴ closed forms below
(4e-13, 1e-13) and the recurrence by finite differences (2e-7 relative).


By powers of Z, h₄₀⁽⁰⁾ depends on ξ only and h₄₀⁽⁴⁾ = -(1+sinα)(2+5sinα)/(9 sinα)
on α only:

    ψ₄₀⁽⁰⁾ = (12E² - 11E + 6a₂₁ + 1)/1152 - (72E a₂₁ + E - 30a₂₁ - 2) sinα cosθ/720
    ψ₄₀⁽⁴⁾ = sinα/36 + 7/288 + 2(4cos²α-1)/(135π)            (pure)

Z¹ to Z³ involve V × (dilogarithmic and elliptic parts of ψ₃₀ and ψ₂₀).

### Reductions of ψ₄₀ (`a21_sector.jl`, `e_sector.jl`)

Write L_k⁺ for the pure inverse of Λ² - k(k+4) and P₄ for the projection onto the
k = 4 harmonics.

* **a₂₁ sector.** ψ₃₀ depends on a₂₁ only through -a₂₁ψ₃₁/(ZB), B = (π-2)/(3π). With
  h₄₁ = 24ψ₄₂ - 2Vψ₃₁ + 2Eψ₂₁ this gives exactly ∂ψ₄₀/∂a₂₁ = -ψ₄₁ᵖᵘʳᵉ/(ZB)
  (checked to 1e-12). Reason: ∂Ψ/∂a₂₁ is itself a local solution, R²Y₂₁ + ….
* **E sector.** (H-E)∂_EΨ = Ψ gives (Λ²-32)∂_Eψ₄₀ = 12∂_Eψ₄₁ - 2V∂_Eψ₃₀ - E/3 + 2ψ₂₀.
  Everything is algebraic except 2Zχ, and the resolvent identity
  L₄⁺χ = (L₄⁺h_χ - χ + P₄χ)/20, with (Λ²-12)χ = h_χ = 2ς/(3ξ sin α) - 8(π-2)v/(3π),
  reduces it to χ itself and the two-centre algebraic source h_χ (checked to 3e-13; the
  P₄χ part is again ∝ 1-2sin²α sin²θ).
* **Algebraic two-centre sources** (polynomial/(d±ξ)) are one-dimensional Feynman
  integrals with the same zonal solutions as ψ₄₁⁽²⁾; sources depending on one centre
  only are ODEs in one angle.

What remains is the pure-Z core L₄⁺[-2V × (Clausen part 𝒟 and elliptic panels of ψ₃₀)].
In the Feynman framework each multiplication by V adds one integration parameter, so
this core comes out as a two-dimensional integral unless one integral is done in
closed form, which is what the ψ₃₀ derivation achieved at k = 3.

### ψ₃₀ panels as superpositions of inverse linear forms (`panel_uform.jl`, `core_a0*.jl`)

With atan(a/b)/a = ∫₀¹ b du/(b² + a²u²), each ψ₃₀ panel is
c₀N₀ + c₂N₂ = (2/π)∫₀¹dt∫₀¹du L(t)(c₀ + c₂t²)(-2yt)/D with D = 4y²t² + P u² linear in
x ∈ S³: D = C(1 - q(t,u)·x), |q| < 1 (1 - |q|² = 4t²u²(1+t²)²/C²). For the F_a panels
(b = ±cos α, -2y = √2 ξ) the factor 1/ξ of V₀ cancels the ξ, so
L₄⁺[-2V₀F_a] is a two-dimensional integral of the same elementary zonal solutions as
ψ₄₁⁽²⁾. Checked against the Green's function: 3e-10 and 4e-10 at two points with a
graded 10-point rule (the integrand has log singularities at u → 0 and t → 1).

Numerical details: per-point zonal solutions are made pure with the closed-form
Funk–Hecke projection M(ρ) (`harmonic_proj.jl`, checked to 1e-16); for ρ < 0.3 they come
from the series f = Σ fₙ(ρc)ⁿ → Chebyshev U_m(c), each divided by 4m(m+2) - 32 with m = 2
dropped (agrees with the closed form to 1e-15 at ρ = 0.25-0.29, where both are accurate).

Remaining core pieces: V₁F_r with matching centre is also two-dimensional; V₀F_r and
V₁F_a need one or two more Feynman parameters (three- and four-dimensional), and the
Clausen part 𝒟 has no linear-form representation yet.

