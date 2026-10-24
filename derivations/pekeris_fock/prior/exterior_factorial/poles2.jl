include("exterior_lib.jl")
using Printf
b = cuspmodel(2); E = cuspsolve(b).E
Ys = riccati_coeffs(b, E, 60; prec = 256)
# asymptotic series at Re R = 14 is accurate to ~1e-16 (see test1)
yasym(R) = sum(ComplexF64.(Ys[n+1]) / R^n for n = 0:40)
U = b.U; Dm = b.D; T = b.W + 2E * I; Lam = Diagonal(b.lam)
f(R, Y) = -(Y * Y + (5 / R * I + 2U) * Y - Lam / R^2 + Dm / R + T)
best = []
for y in 0:0.05:12
    R = 14.0 + im * y; Y = yasym(R); h = -0.004; mx = 0.0; at = 0.0
    while real(R) > 0.6
        k1 = f(R, Y); k2 = f(R + h / 2, Y + h / 2 * k1); k3 = f(R + h / 2, Y + h / 2 * k2); k4 = f(R + h, Y + h * k3)
        Y += h / 6 * (k1 + 2k2 + 2k3 + k4); R += h
        nY = norm(Y); nY > mx && (mx = nY; at = real(R))
        (nY > 1e4 || !isfinite(nY)) && break
    end
    push!(best, (y, mx, at))
end
for (y, mx, at) in best
    mx > 30 && @printf("Im R=%.2f  max|Y|=%.2e at Re R=%.3f\n", y, mx, at)
end
