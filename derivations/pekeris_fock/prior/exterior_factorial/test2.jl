include("exterior_lib.jl")
using Printf
b = cuspmodel(2); E = cuspsolve(b).E
N = 400
for prec in (1024, 3072)
    Ys = riccati_coeffs(b, E, N + 2; prec)
    stir = stirling1(N)
    for omega in (0.22, 0.3, 0.45)
        Bk = factorial_coeffs(Ys, omega, N; stir)
        for R in (1.0, 2.0)
            Yr = cuspexterior(b, E, R)
            ps = factorial_partial_sums(Ys[1], Bk, omega, R)
            @printf("prec=%d omega=%.2f R=%.1f: %s\n", prec, omega, R, join([@sprintf("%.0e", norm(p - Yr)) for p in ps[1:25:end]], " "))
        end
    end
end
