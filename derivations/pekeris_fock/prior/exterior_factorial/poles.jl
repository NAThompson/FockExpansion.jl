# Diagnostic: locate complex singularities of the Riccati solution Y(R) (d=2).
include("exterior_lib.jl")
using Printf
b = cuspmodel(2); E = cuspsolve(b).E
Ys = riccati_coeffs(b, E, 160; prec = 1024)
omega = 0.3; Bk = factorial_coeffs(Ys, omega, 150)
function ycomplex(R)
    setprecision(BigFloat, 1024) do
        x = Complex{BigFloat}(R) / omega; p = x; S = Complex{BigFloat}.(Ys[1])
        for k = 0:length(Bk)-1
            S += Bk[k+1] / p; p *= x + k + 1
        end
        ComplexF64.(S)
    end
end
U = b.U; Dm = b.D; T = b.W + 2E * I; Lam = Diagonal(b.lam)
f(R, Y) = -(Y * Y + (5 / R * I + 2U) * Y - Lam / R^2 + Dm / R + T)
for y in -4:0.5:4
    R = 6.0 + im * y
    Y = ycomplex(R)
    h = -0.002; out = Float64[]
    while real(R) > 0.3
        k1 = f(R, Y); k2 = f(R + h / 2, Y + h / 2 * k1); k3 = f(R + h / 2, Y + h / 2 * k2); k4 = f(R + h, Y + h * k3)
        Y += h / 6 * (k1 + 2k2 + 2k3 + k4); R += h
        if norm(Y) > 1e3 || !all(isfinite, Y)
            @printf("y=%5.2f: blow-up near Re R = %.3f\n", y, real(R)); break
        end
    end
    real(R) <= 0.3 && @printf("y=%5.2f: regular down to Re R=0.3, |Y|=%.2e\n", y, norm(Y))
end
