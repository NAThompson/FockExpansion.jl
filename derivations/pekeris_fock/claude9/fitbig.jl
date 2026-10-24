# BigFloat test of H1: raw Feynman psi41^(2) (R = 1 slice) in span of polynomials(r1,r2,r12) x order-2 atoms (+P twists) x optional multiplier.
using LinearAlgebra, Printf
using ClausenFunctions: cl2
setprecision(BigFloat, 256)
function atoms(r1, r2, r12)
    R2 = r1^2 + r2^2
    P = sqrt(2R2 - r12^2)
    d, s = r1 - r2, r1 + r2
    L5 = r12 * P + s * d
    α = 2atan(r2 / r1)
    β = asin((R2 - r12^2) / R2)
    pp = BigFloat(π)
    DL = cl2(α - β) + cl2(pp + β - α) - cl2(pp - α - β) - cl2(α + β)
    lg = β * log((2r1 * r2)^2 / L5^2) + (β + α) * log((s^2 - r12^2) / (2r1 * r2)) + (β - α) * log((r12^2 - d^2) / (2r1 * r2))
    Ω = lg - DL
    l1, l2, l3 = log(sqrt(R2) / s), log((s + r12) / sqrt(R2)), log((d + r12)^2 / L5)
    base = [one(r1), α, β, l1, l2, l3, Ω]
    vcat(base, P .* base)
end
function monos(r1, r2, r12, n)
    v = typeof(r1)[]
    for j = 0:1, i = 0:n, k = 0:n
        i + j + k <= n && push!(v, r1^i * r2^j * r12^k)
    end
    v
end
function load(file)
    pts = Tuple{BigFloat,BigFloat}[]; vals = BigFloat[]
    for l in eachline(file)
        a, t, v = parse.(BigFloat, split(l))
        push!(pts, (a, t)); push!(vals, v)
    end
    pts, vals
end
xi(a, t) = sqrt(2sin(BigFloat(π) / 4 - a / 2)^2 + 2sin(a) * sin(t / 2)^2)
function runfit(pts, vals, sel, n; mult = (r1, r2, r12) -> one(r1), ntest = 200)
    geo = [(cos(a / 2), sin(a / 2), xi(a, t)) for (a, t) in pts]
    rows = map(geo) do (r1, r2, r12)
        A = atoms(r1, r2, r12)[sel]
        m = monos(r1, r2, r12, n)
        vcat([a .* m for a in A]...)
    end
    M = permutedims(hcat(rows...))
    y = vals .* [mult(g...) for g in geo]
    ntr = length(y) - ntest
    # column scaling, then QR least squares on the training rows; residual on held-out rows
    sc = [norm(M[1:ntr, j]) for j in 1:size(M, 2)]
    Ms = M ./ sc'
    F = qr(Ms[1:ntr, :])
    c = F \ y[1:ntr]
    rtr = norm(Ms[1:ntr, :] * c - y[1:ntr]) / norm(y[1:ntr])
    rte = norm(Ms[ntr+1:end, :] * c - y[ntr+1:end]) / norm(y[ntr+1:end])
    @printf("sel=%s n=%d unknowns=%d  train %.2e  test %.2e\n", string(sel), n, size(M, 2), Float64(rtr), Float64(rte))
    c ./ sc, M, y
end
if abspath(PROGRAM_FILE) == @__FILE__
    pts, vals = load(ARGS[1])
    println("points: ", length(vals))
    for (sel, n) in ((1:1, 8), (1:7, 4), (1:7, 6), (1:7, 8), (1:14, 6))
        runfit(pts, vals, sel, n)
    end
end
