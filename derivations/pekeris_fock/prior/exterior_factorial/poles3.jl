# Track poles of the Moebius charts (Y - c)^{-1} in the right half-plane (d=2).
include("exterior_lib.jl")
using Printf
b = cuspmodel(2); E = cuspsolve(b).E
Ys = riccati_coeffs(b, E, 60; prec = 256)
yasym(R) = sum(ComplexF64.(Ys[n+1]) / R^n for n = 0:40)
U = b.U; Dm = b.D; T = b.W + 2E * I; Lam = Diagonal(b.lam)
f(R, Y) = -(Y * Y + (5 / R * I + 2U) * Y - Lam / R^2 + Dm / R + T)
cs = [-8.0, -4.0, -2.0, -1.0, 1.0, 2.0, 4.0, 8.0]
worst = Dict(c => (0.0, 0.0im) for c in cs)
worstY = (0.0, 0.0im)
for y in 0:0.02:12
    R = 14.0 + im * y; Y = yasym(R); h = -0.004
    while real(R) > 0.25
        k1 = f(R, Y); k2 = f(R + h / 2, Y + h / 2 * k1); k3 = f(R + h / 2, Y + h / 2 * k2); k4 = f(R + h, Y + h * k3)
        Y += h / 6 * (k1 + 2k2 + 2k3 + k4); R += h
        nY = norm(Y)
        (nY > 1e3 || !isfinite(nY)) && break
        nY > worstY[1] && real(R) > 0.5 && (global worstY = (nY, R))
        for c in cs
            w = opnorm(inv(Y - c * I))
            w > worst[c][1] && real(R) > 0.5 && (worst[c] = (w, R))
        end
    end
end
println("Y: ", worstY)
for c in cs
    @printf("c=%5.1f  max|(Y-c)^-1| over Re R>0.5 = %.2e at R=%.3f%+.3fi\n", c, worst[c][1], real(worst[c][2]), imag(worst[c][2]))
end
