# How far out does the Fock series reach?

Morgan (1986) proved the Fock expansion converges pointwise for all ρ. These scripts measure
how many terms that takes in practice for the helium ground state (Z = 2, exact E). They use
the Sylvester chain (`../k4/fock_spectral.jl`).

* `growth.jl`, `resolution.jl`: the largest |ψ_kp| over the angles, with all free constants
  zero, falls from 3e-2 at k = 10 to 4e-9 at k = 24. The ratio of successive terms keeps
  shrinking, so the decay is faster than geometric. Grids with n = 48 to 96 agree to about
  1e-11 relative for every k ≤ 24.
* `hylleraas.jl`: a two-exponent Hylleraas reference, s^l t^{2m} u^n e^{-ζs} with
  ζ = 1.6 and 6.0, l+2m+n ≤ 12, N = 504. E is within 2.8e-12 of the exact value.
* `fit.jl K n ρfit`: Ψ = λ(Φ₀₀ + Σ aⱼΦⱼ) is linear in (λ, λaⱼ). The free constants are
  fitted to the Hylleraas function on ρ ≤ ρfit, at 100 angles × 20 radii.
* `density.jl`: P(ρ ≤ 1, 1.5, 2, 2.5, 3) = 0.26, 0.61, 0.84, 0.94, 0.98.

Fit results: maximum error relative to max|ψ| at each ρ.

| K | ρfit | ρ=0.1 | 0.5 | 1.0 | 1.5 | 2.0 | 2.5 |
|---|---|---|---|---|---|---|---|
| 16 | 1.0 | 1.3e-5 | 8e-7 | 1.7e-6 | 4.9e-4 | 8e-2 | 6 |
| 24 | 1.5 | 1.3e-5 | 9e-7 | 5e-7 | 1.8e-6 | 3e-3 | 1.7 |
| 32 | 1.5 | 1.3e-5 | 9e-7 | 5e-7 | 6.5e-7 | 9e-2 | 3e2 |
| 32 | 2.0 | 1.3e-5 | 9e-7 | 6e-7 | 6.8e-7 | 7.3e-7 | 1.7e-2 |

* Inside the fit region, K = 32 reproduces the Hylleraas function to its own accuracy
  (about 7e-7) out to ρ = 2, a sphere that holds 84% of the probability. At ρ = 0.1 the
  1.3e-5 is the Hylleraas error: its basis has no logarithms.
* Outside the fit region the series fails quickly, and it fails faster for larger K. The
  high-order free constants are fixed only by data at large ρ, which is where they matter.
* The low-order constants are stable across all fits. Spread between fits:

  | constant | value | spread |
  |---|---|---|
  | a₂₁ | 0.4767478 | 6e-8 |
  | a₄₀ | −0.2013783 | 5e-7 |
  | a₄₂ | 0.2265200 | 5e-7 |

  ψ(0) is 4.29460776 × (Hylleraas normalization). a₂₁ is the coefficient of sin α cos θ in
  ψ₂₀; a₄₀ and a₄₂ are the coefficients of Y₄₀ and Y₄₂ in ψ₄₀.
