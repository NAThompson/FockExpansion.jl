# Guess a linear ODE sum_{i<=r} p_i(t) B^(i)(t) = 0, deg p_i <= d, for a far-edge log coefficient B given as a Chebyshev
# series in ty, with t = tan(alpha/4), ty = 16 atan(t)/pi - 5 (side +1 far edge, alpha in (pi, 3pi/2)).
using Serialization, LinearAlgebra, Printf
setprecision(BigFloat, parse(Int, get(ENV, "BITS", "256")))
T = BigFloat
P = T(π)
const OFF = parse(Int, get(ENV, "OFF", "5"))
smul(a, b) = [sum(a[i+1] * b[k-i+1] for i = 0:k) for k = 0:length(a)-1]
function sinv(a)
    b = zeros(T, length(a)); b[1] = 1 / a[1]
    for k = 1:length(a)-1
        b[k+1] = -sum(a[i+1] * b[k-i+1] for i = 1:k) / a[1]
    end
    b
end
sint(a, c0) = vcat([c0], [a[k] / k for k = 1:length(a)-1])
# Taylor coefficients of B(t0 + eps) up to order K
function ssqrt(a)   # power-series square root, a[1] > 0
    b = zeros(T, length(a)); b[1] = sqrt(a[1])
    for k = 1:length(a)-1
        b[k+1] = (a[k+1] - sum(b[i+1] * b[k-i+1] for i = 1:k-1; init = zero(T))) / (2b[1])
    end
    b
end
const VAR = get(ENV, "VAR", "t")
function series_at(c, t0, K)
    e = zeros(T, K + 1); e[1] = t0; e[2] = 1
    if VAR == "u"   # u = cos(alpha); far edge side +1: alpha = 2pi - acos(u), ty = 4 alpha/pi - OFF
        d = sinv(ssqrt([1; zeros(T, K)] .- smul(e, e)))     # 1/sqrt(1-u^2) = -d acos/du
        ac = sint(-d, acos(t0))                               # acos(u) series
        ty = (4 / P) .* (-ac); ty[1] += 8 - OFF
    else
    d = sinv(smul(e, e) .+ [1; zeros(T, K)])
    at = sint(d, atan(t0))
    ty = 16 .* at ./ P; ty[1] -= OFF
    end
    b1 = zeros(T, K + 1); b2 = zeros(T, K + 1)
    for k = length(c)-1:-1:1
        b1, b2 = 2 .* smul(ty, b1) .- b2 .+ [c[k+1]; zeros(T, K)], b1
    end
    smul(ty, b1) .- b2 .+ [c[1]; zeros(T, K)]
end
function guess(c; rmax = 10, maxunk = 120, npts = 0, tlo = parse(Float64, get(ENV, "TLO", "1.15")), thi = parse(Float64, get(ENV, "THI", "2.3")), tol = parse(Float64, get(ENV, "TOL", "1e-18")), inhom = false)
    for tot = 2:maxunk, r = 1:rmax
        nb = inhom ? r + 2 : r + 1
        tot % nb == 0 || continue
        d = tot ÷ nb - 1
        np = npts == 0 ? 2tot + 10 : npts
        ts = [T(tlo) + (T(thi) - T(tlo)) * (i - 1) / (np - 1) for i = 1:np]
        rows = Vector{Vector{T}}()
        for t0 in ts
            s = series_at(c, t0, r)
            der = [s[i+1] * factorial(big(i)) for i = 0:r]
            fun = inhom ? vcat(der, [one(T)]) : der
            push!(rows, [(t0 - (T(tlo) + T(thi)) / 2)^j * f for f in fun for j = 0:d])
        end
        A = permutedims(hcat(rows...))
        A ./= [norm(A[:, k]) for k = 1:size(A, 2)]'
        F = qr(A, ColumnNorm())
        dg = abs.(diag(F.R))
        q = Float64(dg[end] / dg[1]); q2 = Float64(dg[end-1] / dg[1])
        if q < tol && q2 > 1e8 * q
            @printf("  ODE r=%d d=%d inhom=%s (%d unknowns): last %.1e prev %.1e\n", r, d, inhom, tot, q, q2)
            return (r, d, q, q2)
        end
    end
    println("  none up to ", maxunk, " unknowns")
    nothing
end
if abspath(PROGRAM_FILE) == @__FILE__
    L = deserialize(get(ENV, "IN", joinpath(@__DIR__, "logcoefs_side1.jls")))
    sel = get(ENV, "SEL", "")
    maxunk = parse(Int, get(ENV, "MAXUNK", "80"))
    for k in sort(collect(keys(L)); by = string)
        isempty(sel) || occursin(sel, string(k)) || continue
        println(k)
        guess(L[k]; maxunk, inhom = get(ENV, "INHOM", "0") == "1")
    end
end
