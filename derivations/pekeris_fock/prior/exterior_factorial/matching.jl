# Replace cuspexterior in the matching by: inverse factorial series for Y at R1,
# then one Taylor expansion about R1 evaluated at rm = 1. Solve for E and a21.
include("exterior_lib.jl")
using Printf

function series_exterior(b, E, rm; R1, K, order, prec = 768)
    slow, fast = pencil_roots(b, E)
    g = minimum(real(fast)) - maximum(real(slow))
    omega = 1.15 * log(2) / g
    Ys = riccati_coeffs(b, E, K + 2; prec)
    Bk = factorial_coeffs(Ys, omega, K)
    Y1 = factorial_Y(Ys, Bk, omega, R1)
    taylor_connect(b, E, Y1, R1, rm, order)[end]
end

function series_mismatch(b, E; rm = 1.0, details = false, kw...)
    N = max(44, maximum(b.K) + 28)
    c = cuspcoeff(b, E, N); F, D, DD = evalseries(c, rm)
    M = D / F - series_exterior(b, E, rm; kw...)
    ev = eigen(M); idx = argmin(real(ev.values)); val = real(ev.values[idx])
    details || return val
    a = real(F \ ev.vectors[:, idx]); a /= a[1]
    a21 = (length(a) >= 2 ? 2a[2] : 0.0) + b.Z * b.g * (π + 4) / (6π) - b.g^2 / 8
    (; E, a21, match = val)
end

function series_solve(b, E0; kw...)
    lo, hi = E0 - 2e-4, E0 + 2e-4
    f = series_mismatch(b, lo; kw...); g = series_mismatch(b, hi; kw...)
    f * g < 0 || error("no bracket $f $g")
    for _ = 1:30
        mid = (lo * g - hi * f) / (g - f)
        v = series_mismatch(b, mid; kw...)
        abs(v) < 2e-11 && return series_mismatch(b, mid; details = true, kw...)
        sign(v) == sign(f) ? (lo = mid; f = v) : (hi = mid; g = v)
    end
    error("no convergence")
end

ds = isempty(ARGS) ? [2, 4, 6] : parse.(Int, ARGS)
for d in ds
    b = cuspmodel(d)
    z = cuspsolve(b)
    settings = d == 2 ? (R1 = 4.0, K = 160, order = 220) : d == 4 ? (R1 = 6.0, K = 260, order = 320) : (R1 = 9.0, K = 260, order = 420)
    t = @elapsed s = series_solve(b, z.E; settings...)
    @printf("d=%d panels:  E=%.13f a21=%.12f\n", d, z.E, z.a21)
    @printf("d=%d series:  E=%.13f a21=%.12f match=%.1e  (|dE|=%.1e |da21|=%.1e, %.0fs, %s)\n",
        d, s.E, s.a21, s.match, abs(s.E - z.E), abs(s.a21 - z.a21), t, string(settings))
    flush(stdout)
end
