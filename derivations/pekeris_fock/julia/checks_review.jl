# Response to Codex review item 3: S asymmetry, eigen-residual, Gram positivity.
using Printf
include("galerkin_log.jl")
for N in (8, 10, 12), v in (:none, :fixedlog, :fixedlog2)
    r = solve(N, v)
    @printf("%-9s N=%2d size=%4d  asymH=%.1e asymS=%.1e  resid=%.1e  S posdef=%s\n", v, N, r.n,
        Float64(r.asym), Float64(r.asymS), Float64(r.resid), r.posdef)
end
