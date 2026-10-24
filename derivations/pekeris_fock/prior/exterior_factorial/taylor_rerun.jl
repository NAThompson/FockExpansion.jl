# Two-series exterior (factorial series at R1, one Taylor expansion about R1) for d=4,6
# with R1 inside the fast-convergence region of the factorial series.
include("exterior_lib.jl")
using Printf
for (d, R1, K, order) in ((4, 6.0, 260, 360), (6, 9.0, 260, 460))
    b = cuspmodel(d); E = cuspsolve(b).E
    slow, fast = pencil_roots(b, E); g = minimum(real(fast)) - maximum(real(slow))
    omega = round(1.15 * log(2) / g, digits = 3)
    Ys = riccati_coeffs(b, E, K + 2; prec = 1024)
    Bk = factorial_coeffs(Ys, omega, K)
    Y1 = factorial_Y(Ys, Bk, omega, R1)
    @printf("d=%d R1=%.0f factorial error at R1: %.1e\n", d, R1, norm(Y1 - cuspexterior(b, E, R1; rmax = 30.0, order = 30, fraction = 0.3)))
    open(joinpath(@__DIR__, "taylor_d$d.csv"), "w") do io
        println(io, "R,order,err")
        for R0 in (0.75, 1.0, 2.0)
            Yr2 = cuspexterior(b, E, R0; rmax = 30.0, order = 30, fraction = 0.3)
            errs = [norm(y - Yr2) for y in taylor_connect(b, E, Y1, R1, R0, order)]
            foreach(n -> println(io, "$R0,$n,$(errs[n])"), eachindex(errs))
            @printf("  R0=%.2f: err(100,200,%d)=%.1e %.1e %.1e first<1e-12 at order %s\n", R0, order, errs[100], errs[200], errs[end], string(findfirst(<(1e-12), errs)))
        end
    end
    flush(stdout)
end
