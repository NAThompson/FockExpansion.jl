# d=2: inverse factorial series of the Plucker vector xi vs of the Y chart, same omega/precision.
include("plucker_lib.jl")
using Printf
const PREC = 3072
const N = 800
b = cuspmodel(2); E = cuspsolve(b).E; m = 4
Rs = (0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0)
refs = Dict(R => cuspexterior(b, E, R; rmax = 30.0, order = 30, fraction = 0.3) for R in Rs)
for R in Rs
    @printf("reference check R=%.2f: |panel default - panel refined| = %.1e\n", R, norm(cuspexterior(b, E, R) - refs[R]))
end
t1 = @elapsed ps = plucker_series(b, E, N + 1; prec = PREC)
t2 = @elapsed Ys = riccati_coeffs(b, E, N + 2; prec = PREC)
@printf("series: plucker %.0fs, riccati %.0fs, sigma=%.8f beta=%.8f\n", t1, t2, Float64(ps.sigma), Float64(ps.beta))
stir = stirling1(N)
"""Least-squares slope of log(err) vs log(K) over K in [400, 800] -> abscissa p = R + omega*slope."""
function abscissa(errs, R, omega)
    ks = 400:N; x = log.(ks); y = log.(errs[ks.+1])
    s = sum((x .- mean(x)) .* (y .- mean(y))) / sum((x .- mean(x)) .^ 2)
    s, R + omega * s
end
mean(v) = sum(v) / length(v)
open(joinpath(@__DIR__, "plucker_conv.csv"), "w") do io
    println(io, "chart,omega,R,K,err")
    for omega in (0.224, 0.3, 0.4)
        t = @elapsed begin
            Bx = factorial_coeffs_vec(ps.xi, omega, N; prec = PREC, stir)
            By = factorial_coeffs(Ys, omega, N; prec = PREC, stir)
        end
        @printf("omega=%.3f (coefficients %.0fs)\n", omega, t)
        for R in Rs
            Yr = refs[R]
            ref = wedge(big.([I; Yr]), ps.S)
            ex, ey, eg = Float64[], Float64[], Float64[]
            setprecision(BigFloat, PREC) do
                x = big(R) / big(omega); poch = x; s = copy(ps.xi[1])
                for k = 0:N
                    s += Bx[k+1] / poch; poch *= x + k + 1
                    push!(ex, Float64(proj_dist(s, ref)))
                    push!(eg, norm(Float64.(graph_from_plucker(s, ps.S, m)) - Yr))
                end
            end
            ey = [norm(p - Yr) for p in factorial_partial_sums(Ys[1], By, omega, R; prec = PREC)]
            for k = 0:N
                println(io, "plucker_angle,$omega,$R,$k,$(ex[k+1])")
                println(io, "plucker_Y,$omega,$R,$k,$(eg[k+1])")
                println(io, "Ychart,$omega,$R,$k,$(ey[k+1])")
            end
            sx, px = abscissa(ex, R, omega); sy, py = abscissa(ey, R, omega)
            @printf("  R=%.2f plucker angle K=100,200,400,800: %.1e %.1e %.1e %.1e (slope %.2f, p=%.2f) | Y from xi: %.1e\n",
                R, ex[101], ex[201], ex[401], ex[801], sx, px, eg[801])
            @printf("         Y-chart        K=100,200,400,800: %.1e %.1e %.1e %.1e (slope %.2f, p=%.2f)\n",
                ey[101], ey[201], ey[401], ey[801], sy, py)
            flush(stdout)
        end
    end
end
