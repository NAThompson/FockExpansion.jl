# Steklov–Poincaré matching for helium: H5Web files

Files are in `/Users/nathompson7/FockExpansion.jl/examples/sciviz/output/` (gitignored).
Angles follow the repo convention: r1/R = cos(α/2), r2/R = sin(α/2), r12/R = sqrt(1 − sin α cos θ).
Coalescence overlays (magenta) mark r2 = 0 (α = 0), r1 = 0 (α = π) and r12 = 0 (α = π/2, θ = 0).

Pipeline: `compute_*.jl` run in the steklov env (Double64 exterior) and write `data/*.jls`;
`write_*.jl` run in the sciviz env (`julia -t auto --project=examples/sciviz ...`) and write the .nc files.

## steklov_convergence.nc (parsed from the study logs)

Every family is a 2D dataset (x axis, `series`). `multi_dim = series` overlays all curves, and the
`series_label` dataset names them (NaN = not run).

- `/by_d/a21`: a21 against the interior channel degree d. Row 1 is Codex's hyperspherical-channel
  sequence (d = 4..24). The other rows are Steklov runs labelled by (R*, L, exterior basis).
  Overlays: the reference 0.476747879 (`reference`, an hline) and the study estimate
  0.47674788 ± 1e-8 (`steklov_estimate`, a band).
  Look for: at equal d, every Steklov curve sits above Codex's. The Steklov runs continue to d = 48.
- `/by_d/a21_minus_reference`: the same data minus the reference, on a symlog scale (threshold 1e-9).
  Look for: the curves cross zero between d = 24 and d = 48 for L ≥ 10. Low L (4, 6) settles about 1.5e-6 and 5e-7 too high.
- `/by_d/lambda_min`: the smallest matching eigenvalue against d (log scale). It is set mostly by L, not d.
- `/increments/abs_increment_per_unit_d`: |a21(d) − a21(d_prev)| / Δd on log scales (d has a log axis).
  `abs_increment_main` shows Codex against Steklov L = 12, B12P6, with the ∝ d^-6.87 guide drawn as a series.
  Look for: parallel straight lines. The interior channel truncation is algebraic in d and does not depend on L.
- `/by_L/a21`, `/by_L/lambda_min`: a21 and λ_min against the trace degree L, for R* ∈ {0.5, 0.75, 1, 1.5, 2}
  and several d. Look for: a21 falls by geometrically shrinking steps (−4.9e-8, −1.8e-8, −6.2e-9),
  and λ_min falls about 2× per step of 2 in L.
- `/exterior_basis/a21_minus_B12P6`: a21(basis) − a21(B12P6) at R* = 1, for d = 16 and d = 24 (symlog).
  B8P2 is off by about 1e-6. Every other basis is within a few 1e-9 for L ≥ 8.
- `/by_Rstar/a21`: a21 against R* at fixed (L, d). The spread is about 2e-8 at L = 12, d = 40.

## steklov_fields.nc (production point R* = 1, L = 10, d = 32, B10P6, E = E0)

Heatmaps on (α, θ); every map also carries `coordinates = "x y z"`, so it can be drawn on the hypersphere surface.
Global attributes: λ_min, λ_2, a21 and the basis sizes.

- `Phi_trace_cusp`: the matched trace in cusp gauge, Φ(R*, Ω) = Σ_{deg ≤ L} t_i Y_i (smooth by construction).
- `Psi_trace`: the physical trace exp(R* U) Φ, with Ψ(0) = 1. The cusps come from exp(R* U).
- `dRPsi_interior`: ∂_R Ψ at R* from the interior Fock series, using all 289 channels.
- `dRPsi_exterior`: the pointwise ∂_R Ψ at R* of the exterior Ritz solution (36 lifts + 602 corrections: 386 Hylleraas, 216 channel).
- `dRPsi_int_minus_ext`: the flux mismatch for the matched seeds (symlog, with coalescence overlays).
  `_relative` divides it by max|∂_R Ψ_int| (domain ±1, so it shares a colour scale with the wrong seed).
  Look for: about 1% at most, concentrated on the electron–nucleus lines α → 0 and α → π, plus a weak spot at the e–e point.
- `Psi_int_minus_ext`: the value jump at R*. It is 4e-16, because both sides share the trace.
- `*_wrong_seed`: the same quantities for the second matching eigenvector (λ_2 ≈ 7.06), scaled to the same
  trace norm. This is what an unsatisfied boundary condition looks like: an O(1) mismatch over the whole sphere.
- `/trace_projection/*`: each flux projected onto the trace functions g_i = exp(R* U) Y_i (36 harmonics):
  - `int_weak`, `ext_weak`: the fluxes the eigenproblem uses.
  - `int_pointwise`, `ext_pointwise`: the same fluxes by quadrature.
  - `mism_pointwise`: the projected pointwise mismatch.
  - `lambda_M_t`: λ M t.
  `abs_flux_matched` overlays all of them on a log scale.
  Look for: int_pointwise = int_weak to 1e-16. ext_pointwise differs from ext_weak by about 1e-7 to 4e-6.
  So the projected pointwise mismatch (≤ 4e-6) is Ritz derivative error, not the matching residual λ M t (≤ 4e-10).
- `/spectrum`: |λ| of M⁻¹(B_int + A_ext) sorted, Re/Im parts and the first 15 seeds for both eigenvectors.
  Look for: the gap from 3e-8 to 7.06.

## steklov_wavefunction.nc

- `Psi(R, α, θ)`, `Phi`, `log10_abs_Psi`: R is the slider (0.05 … 1.0 from the interior series, 1.1 … 4.0 from the exterior Ritz solution).
  `side` is −1 for interior slices and +1 for exterior slices. `Psi` uses a per-frame colour range; `log10_abs_Psi` uses one range for the whole sequence.
- `Psi_along_R`: Ψ along R at three angle pairs on a log scale, with the R* vline.
- `jump_Psi_at_Rstar` (about 1e-16) and `jump_dRPsi_at_Rstar` (the flux jump, as in the fields file), on the 129² grid (`alpha_f`, `theta_f`).
- `Psi_series_continued_minus_exterior` (slider `R_cont`) and `max_abs_series_minus_exterior`: the interior series
  continued past R* minus the Ritz solution. It grows from 3e-4 at R = 1.05 to 1e5 at R = 1.5. The channels above L are
  pinned to zero trace at R* but have nonzero slope, so they grow like (R/R*)^K outside. The interior series is only valid for R ≤ R*.

## steklov_energy.nc (R* = 1, L = 8, d = 32, B10P6; 11 points in E0 ± 1e-5)

- `lambda_min` against E − E0, with a linear fit (series), λ = 0, and markers at E0 and the zero crossing E*.
  `lambda_min_minus_fit` shows the curvature.
- `lambda_2`: flat at about 7.06 (an isolated mode).
- `a21`: a21 from the λ_min eigenvector against E; it is linear.
- `/L10_check`: the study's three-point L = 10 check (energy_check.log).

## steklov_constants.nc

Parallel 1D datasets over `constant`, labelled by `label` and `summary`. The constants are a21 (physical),
the cusp-gauge seeds a_(0,1), a_(2,0), a_(0,2), a_(2,1), a_(0,3), a_(4,0), a_(2,2), a_(0,4), and the physical
⟨ψ40, Y(2,0)⟩ and ⟨ψ40, Y(0,2)⟩.
- `value_best`: L = 14, B10P6, d = 48, plus the algebraic d-tail.
- `uncertainty`: the maximum of the components below.
- `stable_decimals`: floor(−log10 uncertainty).
- Components: `d_tail`, `L_step_12_to_14`, `basis_B12P6_minus_B10P6`, and `Rstar_spread` (known only for a21, a_(2,0), a_(0,2)).
- `value_L*_d*`: the raw values for each configuration.
