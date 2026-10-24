# Convergence of (a) the inverse factorial series for Y and (b) factorial series at R1
# followed by one Taylor expansion about R1, against cuspexterior at the solved E.
include("exterior_lib.jl")
using Printf
const OUT = @__DIR__
ds = isempty(ARGS) ? [2, 4, 6] : parse.(Int, ARGS)
for d in ds
    b = cuspmodel(d)
    E = cuspsolve(b).E
    slow, fast = pencil_roots(b, E)
    g = minimum(real(fast)) - maximum(real(slow))
    omega = round(1.15 * log(2) / g, digits = 3)
    N = d <= 4 ? 300 : 220
    t = @elapsed Ys = riccati_coeffs(b, E, N + 2; prec = 1024)
    Bk = factorial_coeffs(Ys, omega, N)
    @printf("d=%d m=%d E=%.12f gap=%.4f omega=%.3f N=%d riccati %.1fs\n", d, length(b.K), E, g, omega, N, t)
    open(joinpath(OUT, "factorial_d$d.csv"), "w") do io
        println(io, "R,K,err,referr")
        for R in (0.75, 1.0, 2.0, 3.0, 5.0)
            Yr = cuspexterior(b, E, R)
            Yr2 = cuspexterior(b, E, R; rmax = 30.0, order = 30, fraction = 0.3)
            referr = norm(Yr - Yr2)
            ps = factorial_partial_sums(Ys[1], Bk, omega, R)
            for (k, p) in enumerate(ps)
                println(io, "$R,$(k-1),$(norm(p - Yr2)),$referr")
            end
            @printf("  factorial R=%.2f: err(K=50,100,max)=%.1e %.1e %.1e  ref-vs-ref=%.1e\n", R,
                norm(ps[51] - Yr2), norm(ps[101] - Yr2), norm(ps[end] - Yr2), referr)
        end
    end
    R1 = 4.0
    Y1 = factorial_Y(Ys, Bk, omega, R1)
    open(joinpath(OUT, "taylor_d$d.csv"), "w") do io
        println(io, "R,order,err")
        for R0 in (0.75, 1.0, 2.0)
            Yr2 = cuspexterior(b, E, R0; rmax = 30.0, order = 30, fraction = 0.3)
            ys = taylor_connect(b, E, Y1, R1, R0, 260)
            for (n, y) in enumerate(ys)
                println(io, "$R0,$n,$(norm(y - Yr2))")
            end
            errs = [norm(y - Yr2) for y in ys]
            n12 = findfirst(<(1e-12), errs)
            @printf("  factorial(R1=4)+Taylor R0=%.2f: err(60,120,260)=%.1e %.1e %.1e first<1e-12 at order %s\n",
                R0, errs[60], errs[120], errs[end], string(n12))
        end
    end
    flush(stdout)
end
