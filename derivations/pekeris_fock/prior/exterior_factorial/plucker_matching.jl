# d=2 matching at rm with the Plucker exterior: det[ (F;D) | slow frame ] = 0, no Y involved.
include("plucker_lib.jl")
using Printf
const PREC = parse(Int, get(ENV, "PREC", "3072"))
const K = parse(Int, get(ENV, "K", "800"))
const OMEGA = parse(Float64, get(ENV, "OMEGA", "0.224"))
const RM = parse(Float64, get(ENV, "RM", "1.0"))
b = cuspmodel(2); m = 4; S = subsets(2m, m); stir = stirling1(K)

function interior_wedge(E)
    c = cuspcoeff(b, E, RM > 1.5 ? 110 : max(44, maximum(b.K) + 28))
    F, D, _ = evalseries(c, RM)
    w = wedge(big.([F; D]), S)
    w / norm(w)
end
function f_plucker(E)
    setprecision(BigFloat, PREC) do
        ps = plucker_series(b, E, K + 1; prec = PREC)
        Bx = factorial_coeffs_vec(ps.xi, OMEGA, K; prec = PREC, stir)
        x = big(RM) / big(OMEGA); poch = x; s = copy(ps.xi[1])
        for k = 0:K
            s += Bx[k+1] / poch; poch *= x + k + 1
        end
        Float64(plucker_pairing(interior_wedge(E), s / norm(s), S, 2m))
    end
end
function f_panel(E)
    setprecision(BigFloat, 256) do
        Y = cuspexterior(b, E, RM; rmax = 30.0, order = 30, fraction = 0.3)
        xi = wedge(big.([I; Y]), S)
        Float64(plucker_pairing(interior_wedge(E), xi / norm(xi), S, 2m))
    end
end
function secant(f, E0; h = 1e-7, tol = 1e-15)
    a, b_ = E0 - h, E0 + h; fa, fb = f(a), f(b_)
    for _ = 1:15
        c = b_ - fb * (b_ - a) / (fb - fa)
        a, fa = b_, fb; b_, fb = c, f(c)
        abs(b_ - a) < tol && break
    end
    b_, fb
end
z = cuspsolve(b)
@printf("cuspsolve (panels, Y matching, tol 2e-11): E=%.15f\n", z.E)
Ep, fp = secant(f_panel, z.E)
@printf("panel exterior, Plucker pairing:           E=%.15f  f=%.1e\n", Ep, fp)
t = @elapsed (Ex, fx) = secant(f_plucker, z.E)
@printf("Plucker factorial exterior (K=%d, omega=%.3f, %d bits, rm=%.2f): E=%.15f  f=%.1e  (%.0fs)\n", K, OMEGA, PREC, RM, Ex, fx, t)
@printf("|E_plucker - E_panel| = %.1e   |E_plucker - E_cuspsolve| = %.1e\n", abs(Ex - Ep), abs(Ex - z.E))
