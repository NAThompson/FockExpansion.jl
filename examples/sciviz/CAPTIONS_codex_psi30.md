# Captions: two Codex representations of ψ₃₀

Regenerate (in order):

```sh
cd ~/Documents/Codex/2026-09-05/how-much-compu && \
  PYTHONPATH=work/python-deps python3 ~/FockExpansion.jl/examples/sciviz/codex_series_data.py   # ~7 min, CSVs in output/codex_series/
cd ~/FockExpansion.jl
julia -t 4 --project=examples/sciviz examples/sciviz/codex_green.jl    # Green figures + output/codex_green/*.csv
julia -t 4 --project=examples/sciviz examples/sciviz/codex_series.jl   # series figures, comparison, output/psi30_codex.nc
```

All errors are against the package's elliptic `psi30(α, θ; rtol = 1e-13)` with the He ground-state parameters in `common.jl`.

| File (output/) | Caption | Script |
|---|---|---|
| `codex_green_kernel` | Azimuthally integrated Green kernel G(α,θ; α′,θ′) of (Λ² − 15) for the target (0.7, 0.9), as a function of the source point (clipped at ±12; black: G = 0). Right: radial cuts in four directions are straight lines on a log-r axis, so the only singularity is logarithmic, at the target. | `codex_green.jl` |
| `codex_green_source` | The source 10ψ₃₁ − 2Vψ₂₀ + 2Eψ₁₀ without the S³ measure, on a symmetric log scale. It blows up as 1/r at the nuclei (α′ = 0, π, positive) and changes sign to −∞ at the electron–electron coalescence (π/2, 0); black is the zero contour. | `codex_green.jl` |
| `codex_green_integrand` | What `psi30_green` integrates: source × (G + mirror G) × measure on the folded panel α′ ≤ π/2, for a generic target, one near r₁₂ = 0 and one near a nucleus. Red dashed lines are where `method = :rectangular` splits the domain (through the log singularity); white lines are sign changes. | `codex_green.jl` |
| `codex_green_error_maps` | Error of `psi30_green` at rtol = 1e-8 on a 41 × 41 grid for `:rectangular` and `:duffy`, and the wall time per point (4 threads, so times are noisy). Rectangular is typically 1e-12 to 1e-13 and about 10 ms; Duffy is 1e-10 to 1e-11 and about 100 ms. | `codex_green.jl` |
| `codex_green_rtol` | Green-integral error and cost against rtol at six probe points (circles/solid: rectangular, triangles/dashed: Duffy). The error stays 2–4 decades below rtol and reaches 1e-15 to 1e-14 at rtol = 1e-12, so the two representations agree to round-off. | `codex_green.jl` |
| `codex_series_assembly` | How the Codex analytic series assembles ψ₃₀: elementary pieces, incomplete-beta partial waves (F₂b, F₂c, F₂e), Gauss power series in the nuclear (d₂) and electron–electron (d₁) frames, and point-independent constants C_l from closed forms (l ≤ 1) and the certified J-multipole recurrences (l ≥ 2). | `codex_series.jl` |
| `codex_series_terms` | Size of the l-th Legendre term of each series at four points. Each decays like tan(x/2)ˡ (dotted), with x the polar angle in that series' frame: nuclear for F₂c, F₂e and d₂, electron–electron for d₁. Where tan(x/2) → 1 the decay is algebraic. | `codex_series.jl` |
| `codex_series_convergence` | Left: error of the series truncated at L at seven probe points; dotted lines are the predicted ρᴸ, ρ = max(tan(x/2), tan(y/2)). Right: on the line α = π/2, where the nuclear-frame rate is 1, the monotone envelope of the error falls like L⁻² out to L = 200, and only at L ≳ 1/θ near the coalescence. | `codex_series.jl` |
| `codex_series_inner_N` | Error from truncating the inner Gauss series at N (with L = 40), against N = 300. The error falls like 2^(−N/2), because w = sin(x/2) ≤ 1/√2, and reaches round-off by N ≈ 100 at every probe. | `codex_series.jl` |
| `codex_series_error_maps` | Error of the series at L = 2, 5, 10, 20, 40 on the 41 × 41 grid, beside the predicted log₁₀ ρ⁴¹. The slow cross (α = π/2, and α → 0, π at θ ≈ π/2, where the e-e frame rate tends to 1) matches the prediction. | `codex_series.jl` |
| `codex_r1_recurrence` | The certified complete-integral sequence R₁(l) built exactly from R₁(2) = 29/192, R₁(3) = 1/24 and the step-two first difference d_l ∝ l^(−3/2). Numerically R₁(l) + log 2 ∝ l^(−1/2) → 0; a fit gives R₁(∞) = −log 2 to 4e-7. An independent 25-digit quadrature agrees in every printed digit. | `codex_series.jl` |
| `codex_psi30_comparison` | Accuracy against cost for the three evaluators: median error over six probe points, with min–max bars. Elliptic: 1.7 µs. Green: 2 ms to 0.4 s, reaching 1e-15. Codex series (mpmath): 10–150 ms per point once the constants are cached (3 s for one cold call); it is limited by slow partial-wave convergence near α = π/2, not by cost. | `codex_series.jl` |

NetCDF `psi30_codex.nc` contains the reference ψ₃₀, the Green error and time maps for both methods, the series error at L = 10 and L = 40, the predicted rate ρ, and the stack `log10_err_series(L, theta, alpha)`, which the viewer shows with an L slider.
