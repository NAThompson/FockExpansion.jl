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
ψ⁽²⁾ has two-centre (±e₁, e₂) sources, the same structure that gives ψ₃₀ its
elliptic panels; a closed form for it is open.

## Residual of the truncated series (`residual3.jl`, `residual_study.jl`)

The coefficient of R^(k-2) log^p R in (H-E)Ψ is ½[(Λ²-k(k+4))ψ_kp - h_kp],
with missing coefficients set to zero, so the residual of a truncation is
finite and exact. It reproduces the 256-bit automatic-differentiation values in
`paper/fock_30_coefficient.tex` to 12 digits.

    Ψ⁽³⁾:        R²[s₀ + s₁ log R] + O(R³ log R),  s₀ = Vψ₃₀-Eψ₂₀, s₁ = Vψ₃₁-Eψ₂₁
    +ψ₄₂:        R²[(s₀-ψ₄₂) + (s₁-12ψ₄₂) log R] + …   (removes the Y₄ part of s₁)
    +ψ₄₂+ψ₄₁:    R²[s₀ - 6ψ₄₁ - ψ₄₂] + O(R³ log²R)    (no log at order R²)
    +ψ₄₀:        O(R³ log²R)                           (any Y₄ part of ψ₄₀)

## ψ₄₀ (`psi40_start.py`)

By powers of Z, h₄₀⁽⁰⁾ depends on ξ only and h₄₀⁽⁴⁾ = -(1+sinα)(2+5sinα)/(9 sinα)
on α only:

    ψ₄₀⁽⁰⁾ = (12E² - 11E + 6a₂₁ + 1)/1152 - (72E a₂₁ + E - 30a₂₁ - 2) sinα cosθ/720
    ψ₄₀⁽⁴⁾ = sinα/36 + 7/288 + 2(4cos²α-1)/(135π)            (pure)

Z¹ to Z³ involve V × (dilogarithmic and elliptic parts of ψ₃₀ and ψ₂₀).
