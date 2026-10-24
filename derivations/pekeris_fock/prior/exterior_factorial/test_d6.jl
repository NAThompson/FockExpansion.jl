include("exterior_lib.jl")
using Printf
b = cuspmodel(6); E = cuspsolve(b).E
slow, fast = pencil_roots(b, E); g = minimum(real(fast)) - maximum(real(slow))
Ys = riccati_coeffs(b, E, 262; prec = 1024)
for omega in (0.284, 0.4)
    Bk = factorial_coeffs(Ys, omega, 260)
    for R in (4.0, 6.0, 8.0, 10.0)
        Yr = cuspexterior(b, E, R; rmax = 30.0, order = 30, fraction = 0.3)
        ps = factorial_partial_sums(Ys[1], Bk, omega, R)
        @printf("omega=%.3f R=%.0f: %s\n", omega, R, join([@sprintf("%.0e", norm(p - Yr)) for p in ps[1:20:end]], " "))
    end
end
