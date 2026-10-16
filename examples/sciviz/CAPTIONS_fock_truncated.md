# Truncated Fock expansion: figure captions

All figures use the helium ground-state inputs Z = 2, E = −2.903724377…, a₂₁ = 0.476747879,
and the free k = 4 harmonic constants a₄₀ = a₄₂ = 0. Ψ_K = Σ_{k≤K} R^k Σ_p (log R)^p ψ_{k,p}
with ψ₀₀ = 1. "θ → π" shapes stop at θ = π − 0.03 because `psi40` cannot be evaluated
closer to π (see the report); every other special shape is approached to within 1e-3 rad.

Commands (run from the repository root, in this order):

    julia -t 4 --project=examples/sciviz examples/sciviz/fock_truncated_data.jl    # ~2.5 min, cache
    julia -t 4 --project=examples/sciviz examples/sciviz/fock_truncated_lines.jl   # ~2 min
    julia -t 4 --project=examples/sciviz examples/sciviz/fock_truncated_maps.jl    # ~1 min, needs the cache

| File | Caption | Command |
|---|---|---|
| `output/fock_truncated_along_R.{pdf,png}` | Truncated Fock series Ψ_K(R), K = 0…4, along the hyperradius for four fixed electron-triangle shapes. The truncations agree for R ≲ 0.2 and separate by R ≈ 0.5. Odd truncations turn down and even ones turn up. | `fock_truncated_lines.jl` |
| `output/fock_truncated_terms.{pdf,png}` | Size of each term \|R^k (log R)^p ψ_{k,p}\| in the series at four shapes (solid p = 0, dashed p = 1, dotted p = 2). Near the coalescences and the collinear line, the logarithmic terms ψ₂₁, ψ₃₁ and ψ₄₁ are as large as the non-logarithmic ones. The dips at R = 1 are where log R = 0. | `fock_truncated_lines.jl` |
| `output/fock_truncated_increments.{pdf,png}` | Increments \|Ψ_K − Ψ_{K−1}\| = \|R^K Σ_p log^p R ψ_{K,p}\| on log-log axes. The dashed line marks R_×, where the fourth-order increment first exceeds the third: 1.89 for r₁ = r₂ with θ = π/2, 2.45 near the collinear line, 2.04 near r₁₂ = 0, and only 0.67 near r₂ = 0. | `fock_truncated_lines.jl` |
| `output/fock_truncated_residual.{pdf,png}` | Local-energy error (H − E)Ψ_K/Ψ_K of the truncated series, computed exactly from the recurrence as R^{K−1}Σ_p(Vψ_{K,p} − Eψ_{K−1,p})log^p R − E R^K Σ_p ψ_{K,p} log^p R. It scales as R^{−1}, R⁰, R log R, R² log R and R³ log² R for K = 0…4 (dotted guides). Black markers are an independent 6-D Cartesian Laplacian of Ψ_K computed with ForwardDiff (agreement better than 2e−8 relative). Spikes are zeros of Ψ_K. | `fock_truncated_lines.jl` |
| `output/fock_truncated_cusps.{pdf,png}` | Kato cusp conditions for the truncated series. The direction-averaged log-derivative at the coalescence should equal −Z (electron–nucleus, r₁ → 0) or +½ (electron–electron, r₁₂ → 0). Its error falls as d^K, where d is r₂ or the distance of the pair's midpoint from the nucleus. | `fock_truncated_lines.jl` |
| `output/fock_truncated_r1r2.{pdf,png}` | Ψ_K in the (r₁, r₂) plane at fixed inter-electron angle θ = 0.03 (top: electrons nearly parallel, so the crease along r₁ = r₂ is the electron–electron cusp) and θ = π − 0.03 (bottom). Red: the nodal line Ψ_K = 0, which the truncation introduces. White: contours every 0.1. | `fock_truncated_maps.jl` |
| `output/fock_truncated_r1r2_residual.{pdf,png}` | log₁₀\|(H − E)Ψ_K/Ψ_K\| in the same (r₁, r₂) planes. The region where the local energy is accurate grows from the nucleus with each order. The bright wedge along r₁ = r₂ at θ = 0.03 is the near-coalescence 1/r₁₂, which the low orders cancel only partly. Thin dark curves are sign changes. | `fock_truncated_maps.jl` |
| `output/fock_truncated_shell_residual.{pdf,png}` | log₁₀\|(H − E)Ψ_K/Ψ_K\| over the hyperangles (α, θ) on the shells R ≈ 0.09 and 0.49. The error drops by about an order of magnitude per order at R ≈ 0.09. It stays concentrated at the electron–nucleus edges (α → 0, π) and the electron–electron point (α, θ) = (π/2, 0). | `fock_truncated_maps.jl` |
| `output/fock_truncated_hypersphere.{pdf,png}` | Ψ_K painted on the sphere of electron-triangle shapes, (x, y, z) = (r₁, r₂cos θ, r₂sin θ)/R. Top: build-up K = 1…4 at R ≈ 0.49. Bottom: Ψ₄ normalised to its maximum on the shell R ≈ 0.26, 0.49, 0.75, 1.01. Red dot: r₂ = 0. White dot: r₁₂ = 0. The arc x = 0 is r₁ = 0. | `fock_truncated_maps.jl` |
| `output/fock_truncated_crossing_radius.{pdf,png}` | The radius R_× at which one order's increment first exceeds the previous one, as a function of the shape (α, θ). For 3 vs 2, R_× ranges from 1.15 to 1.86. For 4 vs 3, it ranges from 0.67 at the electron–nucleus edges to 2.45 near the collinear line. This is an empirical guide to where the five-term series stops improving, not a radius of convergence. | `fock_truncated_maps.jl` |

## NetCDF files for H5Web

* `output/fock_truncated.nc` (`fock_truncated_data.jl`):
  * ψ_{k,p}(α, θ) for k ≤ 4 on a 64² grid, with `x y z` hypersphere coordinates.
  * `Psi(K, R, α, θ)`, a 5 × 64³ array with a K slider. The volume view shows the build-up.
  * `local_energy_error(K, R, α, θ)`, the exact (H − E)Ψ_K/Ψ_K.
* `output/fock_truncated_maps.nc` (`fock_truncated_maps.jl`):
  * Ψ_K and the local-energy error on the polar (R, α) grids at θ = 0.03 and θ = π − 0.03, with the matching r₁, r₂ arrays.
  * The crossing-radius maps.
