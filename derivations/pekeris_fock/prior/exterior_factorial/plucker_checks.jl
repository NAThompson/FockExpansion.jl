# Sanity checks: compound derivative, graph recovery, pairing, spectrum gap, formal series.
include("plucker_lib.jl")
using Printf
setprecision(BigFloat, 256)
m = 4; S = subsets(2m, m)
A = big.(randn(2m, 2m)); X = big.(randn(2m, m)); ep = big(1e-30)
C = additive_compound(A, S)
println("compound check: ", Float64(norm((wedge(X + ep * A * X, S) - wedge(X, S)) / ep - C * wedge(X, S))))
Y = big.(randn(m, m))
println("graph check: ", Float64(norm(graph_from_plucker(wedge([I; Y], S), S, m) - Y)))
P = big.(randn(2m, m))
println("pairing check: ", Float64(plucker_pairing(wedge(P, S), wedge(X, S), S, 2m) - det([P X])))
b = cuspmodel(2); E = cuspsolve(b).E
ps = plucker_series(b, E, 60; prec = 512)
ev = eigvals(Float64.(ps.C0))
d = sort(abs.(ev .- Float64(ps.sigma)))
@printf("sigma=%.10f beta=%.10f  |C0 xi0 - sigma xi0|=%.1e\n", Float64(ps.sigma), Float64(ps.beta), Float64(ps.resid0))
@printf("distances from sigma in spec(A0^[4]): %.4f %.4f %.4f %.4f\n", d[1], d[2], d[3], d[4])
slow, fast = pencil_roots(b, E)
@printf("sum slow roots=%.10f, slow-fast gap=%.4f\n", sum(real(slow)), minimum(real(fast)) - maximum(real(slow)))
for R in (8.0, 12.0)
    xs = sum(ps.xi[n+1] / big(R)^n for n = 0:40)
    Yr = cuspexterior(b, E, R; rmax = 30.0, order = 30, fraction = 0.3)
    @printf("R=%.0f asymptotic xi vs panel: proj dist=%.1e  graph err=%.1e\n", R,
        Float64(proj_dist(xs, wedge(big.([I; Yr]), ps.S))), norm(Float64.(graph_from_plucker(xs, ps.S, m)) - Yr))
end
