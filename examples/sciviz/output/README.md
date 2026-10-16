# Output files

The `.nc` files are NetCDF-4 (HDF5) files. Open them in VS Code with the H5Web extension. The PNG and PDF files are the paper figures; the scripts that write them are in `..`.

| File | Question it answers | Open first | Status |
|---|---|---|---|
| `psi30_agreement.nc` | Do the elliptic, Green and Langner forms of ψ₃,₀ agree? (Green to 1e-12; Langner's truncated series is slow near r₁₂ = 0 and θ = π) | root group (panels) | paper figure |
| `psi30.nc` | What does ψ₃,₀ look like, how does it depend on E, Z and a₂₁, and does it satisfy its recurrence? | `psi30` | paper figure |
| `psi30_codex.nc` | How accurate are Codex's Green integral and analytic series for ψ₃,₀? | `log10_err_series` (L slider) | paper figure |
| `fock_truncated.nc` | How does the truncated Fock series Ψ converge in K, and where does it break down in R? | `Psi` (multiples over K) | paper figure |
| `fock_truncated_maps.nc` | At what radius does each higher order overtake the one below it? | `crossing_radius_4_3` | paper figure |
| `hierarchy.nc` | What does the whole hierarchy ψ_{k,p} look like, up to k = 12? | `psi_normalized` (k×p multiples) | paper figure |
| `hierarchy_growth.nc` | How fast do the coefficients grow with k? | `max_abs_psi` | paper figure |
| `hierarchy_spectrum.nc` | How fast do the hyperspherical harmonic amplitudes decay? | `amplitude` | verification |
| `hierarchy_residual.nc` | Does the spectral chain satisfy the recurrence? | `chain_residual` | verification |
| `hierarchy_package.nc` | The package's closed forms for k ≤ 4 (their only copy): root on a 91² grid, `/fock_grid` on the 64² grid of `fock_truncated*` | `psi40` | verification |
| `patched.nc` | Codex's patched global representation: blend, budget, growth bounds | `blend/Z` (with W and F₃) | paper figure |

Run `hierarchy_compute.jl` before `fock_truncated_data.jl` and `fock_truncated_maps.jl`: they read the ψ_{k,p} from `hierarchy_package.nc`.
