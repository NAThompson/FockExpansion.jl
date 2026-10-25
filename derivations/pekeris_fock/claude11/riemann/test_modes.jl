# Check K32(Y;X,eta) = T3(z0)/(2pi sqrt(z0^2-1)) [z0>1] against generic numerical c=32 modes (SOps eigenvectors, n=60).
include(joinpath(@__DIR__, "..", "..", "claude9", "fock_spectral_copy.jl"))
using QuadGK
setprecision(BigFloat, 128); T = BigFloat
n = 60; g = SGrid(T, n); o = SOps(g, 32)
P = Float64(π); b = 3P / 2
xn = [cos(P * (j + 0.5) / n) for j = 0:n-1]; bw = [(-1)^j * sin(P * (j + 0.5) / n) for j = 0:n-1]
function hval(v, z)        # mode vector on Y grid (Y = pi + pi/2 x), value at Y = z
    x = 2(z - P) / P
    ax = bw ./ (x .- xn); dot(ax, v) / sum(ax)
end
z0(Y, X, e) = cos((b + Y) / 2) * cos((X + e) / 2) + cos((Y - b) / 2) * cos((e - X) / 2)
K(Y, X, e) = (z = z0(Y, X, e); z > 1 ? (4z^3 - 3z) / (2P * sqrt(z^2 - 1)) : 0.0)
function lhs(Y, X, v)
    es = range(P / 2, 3P / 2; length = 2001); f(e) = z0(Y, X, e) - 1
    cuts = [P / 2]
    for i = 1:2000
        f(es[i]) * f(es[i+1]) < 0 || continue
        lo, hi = es[i], es[i+1]
        for _ = 1:60; m = (lo + hi) / 2; f(lo) * f(m) <= 0 ? (hi = m) : (lo = m); end
        push!(cuts, (lo + hi) / 2)
    end
    push!(cuts, 3P / 2)
    quadgk(e -> K(Y, X, e) * hval(v, e), cuts...; rtol = 1e-11)[1]
end
λ = Float64.(o.λY); V = Float64.(o.VY)
for j in sortperm(abs.(λ))[[4, 8, 15, 25]]
    v = V[:, j]
    for (X, Y) in ((2.2, 2.9), (3.5, 2.0), (4.3, 4.0))
        r = lhs(Y, X, v) / (hval(v, Y) * hval(v, X) / hval(v, b))
        println("mode λ=", round(λ[j], sigdigits = 6), " (X,Y)=($X,$Y): ratio ", r)
    end
end
