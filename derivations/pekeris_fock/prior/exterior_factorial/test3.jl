include("exterior_lib.jl")
using Printf
b = cuspmodel(2); E = cuspsolve(b).E
N = 300
Ys = riccati_coeffs(b, E, N + 2; prec = 2048)
stir = stirling1(N)
for c in (4.0, 8.0, 16.0)
    Ws = chart_coeffs(Ys, big(c))
    for omega in (0.22, 0.3)
        Bk = factorial_coeffs(Ws, omega, N; stir)
        for R in (0.75, 1.0, 2.0)
            Yr = cuspexterior(b, E, R)
            ps = factorial_partial_sums(Ws[1], Bk, omega, R)
            errs = [norm((c * I - inv(p)) - Yr) for p in ps]
            @printf("c=%4.1f omega=%.2f R=%.2f: %s\n", c, omega, R, join([@sprintf("%.0e", e) for e in errs[1:25:end]], " "))
        end
    end
end
