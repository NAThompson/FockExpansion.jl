# Fock hierarchy figures

All runs: `julia -t 4 --project=examples/sciviz examples/sciviz/<script>.jl`, from the repo root.
Run `hierarchy_compute.jl` first (about 2 min, mostly ψ41 and ψ40). It writes `output/hierarchy.nc`
and `output/hierarchy_package.nc`, which the other scripts read. Shared code is in
`hierarchy_common.jl`. It includes `derivations/k4/fock_spectral.jl` (the spectral chain solver)
without changing it.
Parameters: Z = 2, E = E_He, a₂₁ = 0.47674787900. The free harmonic parts of ψ_k0 for even k ≥ 4
(a₄₀, a₄₂, a₆ₗ, …) are set to zero, as in the package's `psi40`.

| File | Caption | Command |
|---|---|---|
| `output/hierarchy_table.{png,pdf}` | "Periodic table" of the Fock coefficients ψ_kp(α, θ) for k ≤ 6, with rows k and columns p (the power of log R). Each panel is normalised by its own max \|ψ_kp\|, which is printed in its title. The black line is the nodal line ψ_kp = 0. The coefficients for k ≥ 5 come from the spectral chain, which matches the package to about 1e-12 for k ≤ 3. | `hierarchy_table.jl` |
| `output/hierarchy_table_k12.{png,pdf}` | The same table extended to k = 12 (49 coefficients). Structure becomes finer as k grows. The leading-log diagonal ψ_{2m,m} consists of pure hyperspherical harmonics. | `hierarchy_table.jl` |
| `output/hierarchy_graph.{png,pdf}` | The recurrence (Λ² − k(k+4))ψ_kp = 2(k+2)(p+1)ψ_k,p+1 + (p+1)(p+2)ψ_k,p+2 − 2Vψ_k−1,p + 2Eψ_k−2,p drawn as a dependency graph on the (k, p) lattice. The solve order is downward in p within each k. For even (resonant) k, solvability of the p equation fixes the Y_kl part of ψ_k,p+1 (dashed arrows), and the ψ_k0 nodes carry the free harmonic constants. | `hierarchy_graph.jl` |
| `output/hierarchy_growth.{png,pdf}` | (a) max \|ψ_kp\| against k up to k = 32 (n = 80 chain; n = 64 agrees to 2e-11 relative). (b) Ratio and root estimates of the radius of convergence from M_k = Σ_p max\|ψ_kp\|. The ratio grows like 0.59 k^0.58 and does not level off, as expected for a series that converges at every R (Morgan 1986). The rigorous bound of 1/4096 per order from the boundary-matching notes is three to four orders of magnitude more conservative. | `hierarchy_growth.jl` |
| `output/hierarchy_spectrum.{png,pdf}` | (a) Amplitude A_K of each ψ_kp (k ≤ 12) in the singlet S-state harmonics of degree K (Λ² = K(K+4)), relative to its RMS. (b) Relative L² error after truncating the harmonic series at degree K. It decays algebraically, like K^−2.3…−2.4 (A_K ~ K^−3), because of the coalescence singularities. Dark rows are the pure-harmonic leading logs ψ_{2m,m}. The dark cell at K = k in each even-k ψ_k0 row is its zeroed free harmonic part. | `hierarchy_spectrum.jl` |
| `output/hierarchy_residual_package.{png,pdf}` | Recurrence residual \|(Λ²−k(k+4))ψ_kp − h_kp\| of the package evaluators, with Λ² applied by ForwardDiff and scaled by the largest term in the equation. Errors sit at the 1e-16 roundoff level except for a growing edge band in ψ30 (up to 2e-11). There is also one bright column at α = π/2 in ψ20: AD there goes through the removable-limit branch of `chi`, which gives a wrong second derivative (value correct). Right: max residual at three points against BigFloat precision. It follows 2^−bits down to 1e-77 at 256 bits. | `hierarchy_residual.jl` |
| `output/hierarchy_residual_chain.{png,pdf}` | The same residual for every chain coefficient with k ≤ 8 (n = 64 barycentric interpolant, ForwardDiff in (α, θ)). The typical level is 1e-13 to 1e-12, and darker bands follow the nodal lines of ψ. Isolated bright pixels (up to about 1e-6) are points that nearly coincide with a Chebyshev node coordinate, where barycentric derivatives cancel. | `hierarchy_residual.jl` |
| `output/hierarchy_convergence.{png,pdf}` | (a) \|chain − package\| for every coefficient the package provides. It is about 1e-12 for k ≤ 3 and ψ42, 1e-10 to 5e-9 for ψ40 (package step = 0.125), and up to 4e-7 for ψ41 near θ → π and α → 0, π, where the package's ψ41 quadrature loses accuracy. (b, c) Error of the Float64 chain against a 128-bit chain. It converges spectrally to about 1e-13 by n ≈ 24, after which roundoff grows slowly with n. At n = 8 the resonant pairs are not resolved and the solver refuses. | `hierarchy_convergence.jl` |
| `output/hierarchy_logterms.{png,pdf}` | Top: the log coefficients ψ21, ψ31, ψ41 and ψ42 (package). Bottom: ψ41 split into Zψ⁽¹⁾, Z²ψ⁽²⁾ (Green/Feynman quadrature) and Z³ψ⁽³⁾, plus its harmonic part c₀Y₄₀ + c₂Y₄₂, which is fixed by the ψ40 solvability condition (max 0.074, comparable to the other parts). Right: c₀ and c₂ are linear in a₂₁ and E; the dots mark helium (c₀ = 0.0045, c₂ = −0.069). | `hierarchy_logterms.jl` |
| `output/hierarchy_leadinglog.{png,pdf}` | The leading-log tower ψ_{2m,m}, m = 1…6, from the chain. Each is a single combination of degree-2m harmonics (the part outside them is about 1e-12). It changes by less than 4e-14 when E, a₂₁, a₄ₗ and a₆ₗ are changed. Its RMS falls from 0.12 (m = 1) to 6e-14 (m = 6). | `hierarchy_logterms.jl` |

## H5Web files

* `output/hierarchy.nc`: `psi(k, p, theta, alpha)` (k ≤ 12, NaN where p > k/2) with k and p sliders,
  `psi_normalized` (each slice divided by its max \|ψ\|, fixed domain [−1, 1]), and `psi_native`
  (k ≤ 24) on the chain's own Chebyshev grid in the separating coordinates X = α − β, Y = α + β.
  The rectangle's corners there are the coalescences.
* `output/hierarchy_package.nc`: package ψ10 … ψ42 on 91², and ψ40 on 25² (`alpha40`, `theta40`).
* `output/hierarchy_residual.nc`: residual stacks (chain k ≤ 8, package) and the BigFloat study.
* `output/hierarchy_spectrum.nc`: `amplitude(k, p, K)` for line plots against the harmonic degree K.
* `output/hierarchy_growth.nc`: `max_abs_psi(p, k)`, `rms_psi(p, k)` and `M(k)`, for log-scale line plots against k.
