# Figures for ψ₃,₀ (third-order nonlogarithmic Fock coefficient)

Unless stated otherwise: helium ground state, Z = 2, E = −2.903724377034…, a₂₁ = 0.476747879.
Coordinates: r₁/R = cos(α/2), r₂/R = sin(α/2), r₁₂/R = ξ = √(1 − sin α cos θ).
So r₂ = 0 on α = 0, r₁ = 0 on α = π, and r₁₂ = 0 at (α, θ) = (π/2, 0).

Regenerate everything (the grids are cached in `output/psi30.nc` and `output/psi30_agreement.nc`;
delete them to recompute, about 1 min):

    julia -t 4 --project=examples/sciviz examples/sciviz/psi30_compute.jl
    julia -t 4 --project=examples/sciviz examples/sciviz/psi30_figures.jl [name ...]

`name` is one of `map hypersphere parts sectors sensitivity slices cusps zoom agreement
derivatives recurrence quadrature`; with no arguments all figures are drawn.

| File (output/…) | Caption | Command argument |
|---|---|---|
| `psi30_map.pdf` | ψ₃,₀(α, θ) for the helium ground state, with level lines. ψ₃,₀ is negative everywhere, and its minimum, ψ₃,₀ ≈ −3.3658, sits at the electron–electron coalescence (star). Green edges mark the electron–nucleus coalescences r₂ = 0 (α = 0) and r₁ = 0 (α = π). | `map` |
| `psi30_hypersphere.pdf` | Left: ψ₃,₀ painted on the shape hemisphere (r₁/R, (r₂/R)cos θ, (r₂/R)sin θ); the line α = 0 collapses to the pole r₂ = 0 and α = π is the arc r₁ = 0. Right: ψ₃,₀ as a height field over (α, θ), showing the conical cusp at r₁₂ = 0. | `hypersphere` |
| `psi30_parts.pdf` | The three pieces of Eq. (answer): (a) the dilogarithmic term 𝒟, (b) the elementary E, a₂₁ term, (c) the sum of the four elliptic contributions F. (d) log₁₀(Σ\|parts\|/\|ψ₃,₀\|): the parts cancel by at most 0.45 digits, so the decomposition is numerically benign everywhere. | `parts` |
| `psi30_Z_sectors.pdf` | ψ₃,₀ at E = a₂₁ = 0 is a cubic polynomial in Z, Σₙ Zⁿ Aₙ(α, θ); the panels show A₀ … A₃ (obtained from Z = 0, 1, 2, 3 and checked at Z = 4 to 4×10⁻¹⁴). The Z³ coefficient is independent of θ. | `sectors` |
| `psi30_Z_shares.pdf` | Share Zⁿ‖Aₙ‖/Σₘ Zᵐ‖Aₘ‖ of each Z sector in the hyperspherical L² norm (weight sin²α sin θ). The Z³ sector dominates from Z ≈ 1.1 onwards; for helium it carries 60% of the norm. Numerically ‖A₀‖ = 1/72 to 16 digits. | `sectors` |
| `psi30_sensitivity.pdf` | ForwardDiff sensitivities of ψ₃,₀ at the helium parameters: (a) ∂/∂E, (b) ∂/∂a₂₁, (c) ∂/∂Z at fixed E and a₂₁. ψ₃,₀ is affine in E and a₂₁, so (a) and (b) are exact coefficient functions. ∂/∂a₂₁ changes sign near θ = π/2 (it is ∝ sin α cos θ apart from the ξ terms). | `sensitivity` |
| `psi30_slices.pdf` | ψ₃,₀ (solid) and its three parts (dashed) along (a) the collinear line θ → 0, (b) the collinear line θ → π, (c) the line r₁ = r₂. The kink at α = π/2 in (a) is the electron–electron cusp. | `slices` |
| `psi30_cusps.pdf` | Kato cusp conditions at third order. (a) Along α = π/2, ∂ψ₃,₀/∂ξ tends to ψ₂,₀(coalescence)/2 (dashed); the lower panel shows linear convergence down to ξ = 10⁻⁶. (b) The θ-average of ψ₃,₀ satisfies ∂⟨ψ₃,₀⟩/∂s → −Zψ₂,₀(r₂ = 0) with s = r₂/R, again with linear convergence. | `cusps` |
| `psi30_coalescence_zoom.pdf` | Neighbourhood of r₁₂ = 0. (a) ψ₃,₀ minus its coalescence value has circular, cone-like level sets. (b) Subtracting the Kato cone ½ψ₂,₀(coal.) ξ leaves a remainder of size 10⁻² with a nodal V. (c) \|remainder\|/ξ² stays bounded, so the remainder is O(ξ²). | `zoom` |
| `psi30_agreement.pdf` | Absolute difference between the compact evaluator and (a) the Liverts Green-function integral (rtol 10⁻¹⁰), (b) Langner's double series truncated at N = 60, J = 320, on a 65² grid. The Green representation agrees to 10⁻¹⁴–10⁻¹² everywhere. The Langner series agrees to round-off only in the interior and degrades to 4×10⁻² near the collinear and coalescence boundaries, where it converges slowly. | `agreement` |
| `psi30_derivatives.pdf` | (a) Gradient norm \|∇ψ₃,₀\| = 2(ψ_α² + ψ_θ²/sin²α)^{1/2} on the hypersphere; it stays finite at r₁₂ = 0, where only its direction jumps. (b) Λ²ψ₃,₀ from ForwardDiff Hessians (pseudo-log colour). It diverges like 1/ξ at r₁₂ = 0 and like 1/r at the nuclear edges, matching −2Vψ₂,₀. | `derivatives` |
| `psi30_recurrence.pdf` | Residual of (Λ² − 21)ψ₃,₀ = 10ψ₃,₁ − 2Vψ₂,₀ + 2Eψ₁,₀, with the left side differentiated by ForwardDiff on a 513² grid: (a) absolute, (b) relative to max(\|Λ²ψ\|, \|RHS\|). The interior sits at 10⁻¹⁵–10⁻¹¹. The sharp-edged patches mirror the switching between fixed Gauss rules in the elliptic quadrature. Values rise above 10⁻⁶ (white) only within 10⁻³ of the boundaries. | `recurrence` |
| `psi30_quadrature_error.pdf` | Quadrature error estimate returned by `psi30_parts` (rtol = 10⁻¹¹); it is uniform to within a factor 3 across the square. | `quadrature` |

## NetCDF files for H5Web

- `output/psi30.nc` (513² interior grid, ε = 10⁻³; `alpha`, `theta` dimension scales, `x y z`
  hypersphere coordinates): `psi30`, `psi30_classical`, `psi30_state`, `psi30_elliptic`,
  `psi30_quadrature_error`, `log10_cancellation`, `dpsi30_dE`, `dpsi30_da21`, `dpsi30_dZ`,
  `sector_Z0` … `sector_Z3`, `grad_norm`, `lambda2_psi30`, `recurrence_rhs`,
  `recurrence_residual`, `recurrence_residual_abs`, `recurrence_residual_rel`. Each field has
  typed `H5WEB_*` view attributes (symlog with an automatic diverging map for signed
  fields, log for errors) and the coalescence overlays.
- `output/psi30_agreement.nc` (65² grid, ε = 0.02): `psi30`, `psi30_green`, `psi30_langner`,
  `abs_diff_green`, `abs_diff_langner`.
