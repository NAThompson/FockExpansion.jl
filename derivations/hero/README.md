# Hero run: the Fock expansion of helium to high order

Uses the Sylvester chain (`../k4/fock_spectral.jl`) in Double64 (`DoubleFloats`, about 31
digits; this directory has its own `Project.toml`). The eigendecompositions are refined by
Newton's method, and long runs keep only the two latest orders (`sfock(...; keep = 2)`).

## Stage 1: the bare chain to k = 120 (`growth_hi.jl`, `growth_compare.jl`, `reach.jl`)

Z = 2, exact E, free constants zero. Grids of n = 128, 160 and 192 points (192 took about
an hour). `stage1_n192.csv` has max|ψ_k| for each k, the n = 160 vs 192 difference, and
the timing.

* Resolution: n = 160 and 192 agree to about 1e-26 relative at every k ≤ 120. n = 128
  holds 1e-26 up to k ≈ 90.
* Decay: max|ψ_k| is 3e-2 at k = 10, 4e-36 at k = 60 and 1.5e-93 at k = 120. It goes
  roughly like c^k/(k/2)!, faster than geometric.
* Reach (`reach.jl`): the largest term max|ψ_kp| ρ^k |log ρ|^p, and the order beyond
  which all terms are below 1e-20 times that:

  | ρ | largest term | at k | terms < 1e-20·largest beyond |
  |---|---|---|---|
  | 2 | 3e1 | 6 | 63 |
  | 3 | 3e3 | 16 | 84 |
  | 4 | 1e6 | 26 | 106 |
  | 5 | 1e9 | 36 | > 120 |
  | 6 | 2e12 | 48 | > 120 |

  The largest term sets the cancellation. Double64 therefore handles ρ ≲ 6–7; beyond that
  needs more digits and K > 120.

## Stage 2: energy from the decay condition (`dirichlet.jl`, prototype)

Ψ = 0 on the hypersphere ρ = ρ_out. Ψ is linear in the free constants at fixed E, so the
confined-helium energy E_D(ρ_out) is where the column-normalised collocation matrix
[Φⱼ(ρ_out, Ωᵢ)] becomes singular. E_D(ρ_out) decreases to the free-atom energy as ρ_out
grows.
