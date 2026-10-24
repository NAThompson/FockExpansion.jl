include("exterior_lib.jl")
using Printf
d = 2
b = cuspmodel(d)
z = cuspsolve(b)
E = z.E
slow, fast = pencil_roots(b, E)
g = minimum(real(fast)) - maximum(real(slow))
println("gap g = $g, slow=$(real(slow))")
N = 200
t = @elapsed Ys = riccati_coeffs(b, E, N; prec = 1024)
println("riccati coeffs $t s")
Y0f, Y1f, _ = cusptail(b, E)
println("Y0 diff ", norm(Float64.(Ys[1]) - Y0f), " Y1 diff ", norm(Float64.(Ys[2]) - Y1f))
for n in (10, 20, 50, 100, 150, 199)
    @printf("n=%d  ||Y_{n+1}||/(n||Y_n||) = %.5f  (1/g = %.5f)\n", n, Float64(norm(Ys[n+2]) / (n * norm(Ys[n+1]))), 1 / g)
end
Yref = cuspexterior(b, E, 12.0)
ap = asymptotic_partial_sums(Ys, 12.0)
println("asymptotic at R=12: ", [@sprintf("%.1e", norm(a - Yref)) for a in ap[1:5:80]])
for omega in (0.3, 0.5, 1.0)
    K = 150
    Bk = factorial_coeffs(Ys, omega, K)
    for R in (1.0, 2.0, 5.0)
        Yr = cuspexterior(b, E, R)
        ps = factorial_partial_sums(Ys[1], Bk, omega, R)
        @printf("omega=%.2f R=%.1f: %s\n", omega, R, join([@sprintf("%.1e", norm(p - Yr)) for p in ps[1:10:end]], " "))
    end
end
