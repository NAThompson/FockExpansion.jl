# Two convergent series, one determinant: Fock interior and Y-chart inverse factorial
# exterior, both in BigFloat, matched at R* by det(D - Y F) = 0. No panels, no Taylor.
include("exterior_lib.jl")
using Printf

# Generic-typed copy of cuspcoeff (same recurrence as cusp_matching.jl).
function cuspcoeff_big(b, E::BigFloat, N)
    m = length(b.K); pmax = N ÷ 2 + 1
    U = big.(b.U); Dm = big.(b.D); T = big.(b.W) + 2E * I; lam = big.(b.lam)
    c = zeros(BigFloat, N + 1, pmax + 2, m, m)
    for n = 0:N
        src = zeros(BigFloat, pmax + 2, m, m)
        for k = 0:n÷2
            if n >= 1
                src[k+1, :, :] -= (2(n - 1) * U + Dm) * c[n, k+1, :, :] + 2(k + 1) * U * c[n, k+2, :, :]
            end
            n >= 2 && (src[k+1, :, :] -= T * c[n-1, k+1, :, :])
        end
        for i = 1:m
            den = n * (n + 4) - lam[i]
            if den == 0
                for k = n÷2:-1:0
                    c[n+1, k+2, i, :] = (src[k+1, i, :] - (k + 1) * (k + 2) * c[n+1, k+3, i, :]) / ((2n + 4) * (k + 1))
                end
                c[n+1, 1, i, i] = 1
            else
                for k = n÷2:-1:0
                    c[n+1, k+1, i, :] = (src[k+1, i, :] - (2n + 4) * (k + 1) * c[n+1, k+2, i, :] - (k + 1) * (k + 2) * c[n+1, k+3, i, :]) / den
                end
            end
        end
    end
    c
end

"""Fock values F and derivatives D at R (BigFloat)."""
function fock_FD(c, R::BigFloat)
    N = size(c, 1) - 1; m = size(c, 3)
    F = zeros(BigFloat, m, m); D = zeros(BigFloat, m, m); l = log(R)
    for n = 0:N, k = 0:n÷2
        A = c[n+1, k+1, :, :]
        p = l^k; p1 = k > 0 ? k * l^(k - 1) : big(0)
        F += R^n * p * A
        D += R^(n - 1) * (n * p + p1) * A
    end
    F, D
end

"""Y at R from the inverse factorial series, kept in BigFloat."""
function yfactorial_big(b, E::BigFloat, R; K, omega, rprec)
    Ys = riccati_coeffs_bigE(b, E, K + 2; prec = rprec)
    setprecision(BigFloat, rprec) do
        Bk = factorial_coeffs(Ys, omega, K)
        x = big(R) / big(omega); poch = x; S = copy(Ys[1])
        for k = 0:K
            S += Bk[k+1] / poch; poch *= x + k + 1
        end
        S
    end
end

# riccati_coeffs with a BigFloat energy (exterior_lib's version takes Float64 E).
function riccati_coeffs_bigE(b, E::BigFloat, N; prec)
    setprecision(BigFloat, prec) do
        Y, _, _ = cusptail(b, Float64(E))
        U = big.(b.U); T = big.(b.W) + 2E * I; Dm = big.(b.D); Lam = Diagonal(big.(b.lam))
        Y = big.(Y); m = size(Y, 1)
        for _ = 1:12
            Q = Y * Y + 2U * Y + T
            norm(Q) < big(2.0)^(-prec + 20) && break
            Y -= reshape(sylvester_operator(Y + 2U, Y) \ vec(Q), m, m)
        end
        S = lu(sylvester_operator(Y + 2U, Y))
        Ys = [Y, reshape(S \ vec(-(5Y + Dm)), m, m)]
        for n = 2:N
            C = (n - 6) * Ys[n]
            for j = 1:n-1
                C -= Ys[j+1] * Ys[n-j+1]
            end
            n == 2 && (C += Lam)
            push!(Ys, reshape(S \ vec(C), m, m))
        end
        Ys
    end
end

"""Matching determinant and amplitudes at energy E."""
function match(b, E::BigFloat; Rs, K, omega, nfock, prec, rprec)
    setprecision(BigFloat, prec) do
        c = cuspcoeff_big(b, E, nfock)
        F, D = fock_FD(c, big(Rs))
        Y = setprecision(BigFloat, prec) do
            BigFloat.(yfactorial_big(b, E, Rs; K, omega, rprec))
        end
        M = D - Y * F
        f = det(M) / prod(norm(M[:, j]) for j in axes(M, 2))
        # null vector with a[1] = 1 (least squares over all rows)
        a = vcat(big(1), -(M[:, 2:end] \ M[:, 1]))
        a21 = (length(a) >= 2 ? 2a[2] : big(0)) + b.Z * b.g * (big(π) + 4) / (6big(π)) - b.g^2 / 8
        (; f, a, a21)
    end
end

function secant_root(fun, E0; h = 1e-9, tol = 1e-24, maxit = 12)
    a = big(E0) - h; bb = big(E0) + h; fa = fun(a); fb = fun(bb)
    for _ = 1:maxit
        cc = bb - fb * (bb - a) / (fb - fa)
        a, fa = bb, fb; bb, fb = cc, fun(cc)
        abs(bb - a) < tol && break
    end
    bb, (fb - fa) / (bb - a)
end

"""Panel root localised by secant on cuspmismatch (Float64)."""
function panel_root(b, E0)
    f(E) = cuspmismatch(b, E)
    a, bb = E0 - 1e-9, E0 + 1e-9; fa, fb = f(a), f(bb)
    for _ = 1:20
        cc = bb - fb * (bb - a) / (fb - fa)
        a, fa = bb, fb; bb, fb = cc, f(cc)
        (abs(bb - a) < 1e-16 || fb == 0) && break
    end
    bb, cuspmismatch(b, bb; details = true).a21
end

base = Dict(
    2 => (Rs = 4.0, K = 200, omega = 0.224, nfock = 110, prec = 384, rprec = 1024),
    4 => (Rs = 6.0, K = 260, omega = 0.256, nfock = 150, prec = 384, rprec = 1024),
    6 => (Rs = 10.0, K = 240, omega = 0.284, nfock = 210, prec = 384, rprec = 1024),
)
variants(p) = [
    "R*+1" => merge(p, (Rs = p.Rs + 1, nfock = p.nfock + 25)),
    "K+60" => merge(p, (K = p.K + 60,)),
    "omega*1.25" => merge(p, (omega = 1.25p.omega,)),
    "fock+30" => merge(p, (nfock = p.nfock + 30,)),
    "prec+128" => merge(p, (prec = p.prec + 128, rprec = p.rprec + 512)),
]

ds = isempty(ARGS) ? [2, 4, 6] : parse.(Int, ARGS)
open(joinpath(@__DIR__, "twoseries.csv"), "a") do io
    for d in ds
        b = cuspmodel(d); z = cuspsolve(b)
        Ep, a21p = panel_root(b, z.E)
        @printf("d=%d panels: cuspsolve E=%.16f  tightened E=%.16f a21=%.14f\n", d, z.E, Ep, a21p)
        p = base[d]
        t = @elapsed (Er, slope) = secant_root(E -> match(b, E; p...).f, Ep)
        r = match(b, Er; p...)
        @printf("d=%d two-series %s\n   E=%.20f a21=%.16f  (%.0fs)\n   |E-E_panel|=%.1e |a21-a21_panel|=%.1e\n",
            d, string(p), Float64(Er), Float64(r.a21), t, Float64(abs(Er - Ep)), Float64(abs(r.a21 - a21p)))
        println(io, "$d,base,$(Float64(Er)),$(Float64(r.a21)),$(Float64(Er-Ep)),$(Float64(r.a21-a21p)),$t")
        flush(stdout)
        for (name, q) in variants(p)
            # one Newton step from the base root with the base slope, then a21 there
            t = @elapsed begin
                fv = match(b, Er; q...).f
                Ev = Er - fv / slope
                rv = match(b, Ev; q...)
            end
            @printf("   variant %-11s dE=%.1e da21=%.1e  residual f=%.1e (%.0fs)\n", name,
                Float64(Ev - Er), Float64(rv.a21 - r.a21), Float64(rv.f), t)
            println(io, "$d,$name,$(Float64(Ev)),$(Float64(rv.a21)),$(Float64(Ev-Er)),$(Float64(rv.a21-r.a21)),$t")
            flush(stdout)
        end
    end
end
