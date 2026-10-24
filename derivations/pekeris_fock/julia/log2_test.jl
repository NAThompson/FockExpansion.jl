# Does a log^2(sigma) sector help? (variant :fixedlog2 = :fixedlog + log^2 σ × polynomials of degree 4..N-2)
using Printf
include("galerkin_log.jl")
E0 = big"-2.903724377034119598311"
@printf("zeta(3) check: %.30f\n", Float64(zeta_em(3)) - 1.2020569031595942)
for N = 6:parse(Int, get(ARGS, 1, "12")), v in (:fixedlog, :fixedlog2)
    r = solve(N, v)
    @printf("%-10s N=%2d size=%4d  E-E0=%.4e  asym=%.1e\n", v, N, r.n, Float64(r.E - E0), Float64(r.asym))
end
