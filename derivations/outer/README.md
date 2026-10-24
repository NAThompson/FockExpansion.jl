# Outer formal solution of helium in the He⁺(1s) + e channel

The Fock series is built at ρ = 0. This is the analogous object at r₂ → ∞, to study whether
a convergent outer expansion with the same properties exists.

## The recurrence (`outer_series.jl`)

r₁ is the inner electron and R = r₂ the outer one. With ψ = Σ_l f_l(r₁, R) P_l(cos θ₁₂) and
1/r₁₂ = Σ_λ r₁^λ/R^{λ+1} P_λ:

    R f_l = e^{-κR} Σ_k R^{σ-k} g_kl(r₁),   g_kl = e^{-Zr₁} × polynomial(r₁)

    L_l g_kl = κ(k-1) g_{k-1,l} + ½[(σ-k+2)(σ-k+1) - l(l+1)] g_{k-2,l}
               - Σ_{λ≥1, l'} (2l+1)(l λ l'; 0 0 0)² r₁^λ g_{k-λ-1,l'}

with L_l = h_l(r₁) - ε₁ₛ the hydrogenic operator shifted by the He⁺ ground-state energy.

* k = 0 gives g₀₀ = e^{-Zr₁} and κ = √(2(-Z²/2 - E)) = 1.3444.
* k = 1 gives σ = (Z-1)/κ = 0.7438.
* For l = 0, L₀ annihilates e^{-Zr₁}. Its coefficient in g_{k-1,0} is fixed by solvability at
  order k, as the harmonic parts are in Fock's recurrence. There are no logarithms, because
  the factor κ(k-1) never vanishes for k ≥ 2.
* L_l is triangular on r^j e^{-Zr}, so every coefficient is exactly polynomial × e^{-Zr₁}.
  The script computes them by downward recursion in 512-bit arithmetic; 40 orders take 7 s.

## Checks

* `residual.jl`: (H-E)ψ_K/ψ_K falls like R^{-(K+1)} at fixed K. At R = 40 it is 2e-13 with
  K = 12. At R = 10 it is smallest near K ≈ 8 and then grows, the usual optimal
  truncation of an asymptotic series.
* `check_hylleraas.jl`: against the Hylleraas reference at r₂ = 4-5, with the amplitude
  fitted at one point, the error falls with K to 2-5e-4. At r₂ ≥ 6 the reference is
  unreliable: its basis decays like e^{-1.6s} and cannot reproduce the e^{-1.344R} tail.

## The series diverges, and why (`growth.jl`, `signs.jl`)

* g_kl grows like Γ(k+β)/A^k with A ≈ κ₂ - κ₁ = 0.846, where κ₂ = √(2(-Z²/8 - E)) is the
  decay rate of the He⁺(n = 2) channel. The ratios g_{k+2}/g_k at k = 36-38 match
  (k+β)(k+β+1)/A² with β ≈ -2.
* The coefficients keep one sign (all negative at small r₁ for k ≤ 38; all positive at
  r₁ = 3). The Borel transform therefore has its singularity on the positive real axis,
  at ζ = κ₂ - κ₁.
* So the 1s series is not Borel summable along the real axis. The two lateral sums differ
  by a multiple of e^{-(κ₂-κ₁)R} relative to the 1s term: exactly the n = 2 channel's own
  decaying solution.

Consequence: no single-channel outer series can converge to the wavefunction. A convergent
outer representation has to be a multi-channel transseries:

    ψ_out = Σ_n C_n e^{-κ_n R} R^{σ_n} (series_n)

Each series is resummed, and the channels are linked by Stokes constants. The channel
amplitudes C_n play the role on the outside that Fock's free constants play inside.
Whether that sum converges is the open question: κ_n → √(-2E) as n → ∞, and the He⁺
continuum lies above that.
