# Is the ladder piece phi3^A (and B, C) in the span of the 115 order-3 atoms x polynomials? 256-bit fit on Codex x19 points.
include(joinpath(@__DIR__, "logcoef.jl"))
using DelimitedFiles
setprecision(BigFloat, 256); T = BigFloat
X19 = joinpath(@__DIR__, "..", "codex", "phase9")
ids = parse.(Int, readlines(joinpath(X19, "x19_fit40_big_ids.txt")))
lines = readlines(joinpath(@__DIR__, "..", "claude9", "spec40.txt"))
pts = [parse.(T, split(lines[i]))[1:2] for i in ids]
A = readdlm(joinpath(X19, "x19_fit40_big_atoms.txt"), BigFloat)
geometry(a, t) = (cos(a / 2), sin(a / 2), sqrt(2sin(T(π) / 4 - a / 2)^2 + 2sin(a) * sin(t / 2)^2))
geo = [geometry(p...) for p in pts]
n = 80; g = SGrid(T, n); P = T(π)
include(joinpath(@__DIR__, "..", "claude9", "psi40_big.jl"))   # interp(g, F, α, θ)
wv0 = 2sqrt(T(2)) .* sin.(g.α) .* cos.(g.β ./ 2); wv1 = 4 .* sin.(g.β) .* (cos.(g.α ./ 2) .+ sin.(g.α ./ 2))
function ladder(lad)
    ψ = Dict{Tuple{Int,Int},Matrix{T}}((0, 0) => ones(T, n, n)); get0(k, p) = get(ψ, (k, p), zeros(T, n, n))
    for k = 1:3
        o = SOps(g, k * (k + 4)); H = iseven(k) ? sharmonics(g, k) : Matrix{T}[]
        wVk = lad[k][1] .* wv0 .- lad[k][2] .* wv1
        for p = (k÷2):-1:0
            wh(f) = g.w .* (2(k + 2) * (p + 1) * f + (p + 1) * (p + 2) * get0(k, p + 2)) - 2wVk .* get0(k - 1, p)
            if iseven(k) && p < k ÷ 2 && haskey(ψ, (k, p + 1))
                G = [sip(g, a, b) for a in H, b in H]; s = [sum(g.W ./ g.w .* h .* wh(get0(k, p + 1))) for h in H]
                ψ[(k, p + 1)] += sum(c * h for (c, h) in zip(-(G \ s) ./ (2(k + 2) * (p + 1)), H))
            end
            ψ[(k, p)] = ssolve(g, o, wh(get0(k, p + 1)), H)
        end
    end
    ψ
end
paths = Dict("A" => [(0, -1), (1, 0), (1, 0)], "B" => [(1, 0), (0, -1), (1, 0)], "C" => [(1, 0), (1, 0), (0, -1)])
Y = zeros(T, length(pts), 4)
for (col, key) in enumerate(("A", "B", "C"))
    ψ = ladder(paths[key])
    Y[:, col] = [interp(g, ψ[(3, 0)], p[1], p[2]) for p in pts]
end
Y[:, 4] = Y[:, 1] + Y[:, 2] + Y[:, 3]
monsb(x, y, c, d) = [x^i * y^j * c^k for j = 0:1 for i = 0:d for k = 0:d if i + j + k <= d]
ntr = 800
for (fam, cols) in (("all115", 1:115), ("order2+classical", 1:23))
    for d in (1, 2, 3)
        nm = length(monsb(geo[1]..., d)); M = zeros(T, length(geo), nm * length(cols))
        for i = 1:length(geo)
            m = monsb(geo[i]..., d)
            for (ai, a) in enumerate(cols), j = 1:nm
                M[i, (ai-1)*nm+j] = A[i, a] * m[j]
            end
        end
        size(M, 2) > ntr - 50 && continue
        for j = 1:size(M, 2); s = norm(M[1:ntr, j]); s > 0 && (M[:, j] ./= s); end
        F = qr(M[1:ntr, :], ColumnNorm()); dd = abs.(diag(F.R)); r = findfirst(<(dd[1] * big"1e-45"), dd); r = isnothing(r) ? length(dd) : r - 1
        yy = F.Q' * Y[1:ntr, :]; co = F.R[1:r, 1:r] \ yy[1:r, :]; res = M[:, F.p[1:r]] * co - Y
        for (j, nmj) in enumerate(("phi3A", "phi3B", "phi3C", "sum=-psi30Z1"))
            te = norm(res[ntr+1:end, j]) / norm(Y[ntr+1:end, j])
            println(fam, " d=", d, " cols=", size(M, 2), " rank=", r, " ", nmj, " held-out rel = ", Float64(te))
        end
        flush(stdout)
    end
end
