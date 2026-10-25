# Physical ladder phi1..phi4 (phi_k = T_k[-2 wV_k phi_{k-1}] with wV_k selected) restricted to the physical edge X = side*pi/2,
# as Chebyshev series in ty; saved for ODE guessing (OFF = 3 for side = +1: alpha = 3pi/4 + pi ty/4, t = tan(alpha/4)).
include(joinpath(@__DIR__, "logcoef.jl"))
setprecision(BigFloat, parse(Int, get(ENV, "BITS", "256"))); T = BigFloat
n = parse(Int, get(ENV, "NGRID", "80")); side = parse(Int, get(ENV, "SIDE", "1"))
lad = eval(Meta.parse(get(ENV, "LADDER", "[(1,0),(1,0),(0,-1),(1,0)]")))
g = SGrid(T, n); P = T(π)
wv0 = 2sqrt(T(2)) .* sin.(g.α) .* cos.(g.β ./ 2); wv1 = 4 .* sin.(g.β) .* (cos.(g.α ./ 2) .+ sin.(g.α ./ 2))
ψ = Dict{Tuple{Int,Int},Matrix{T}}((0, 0) => ones(T, n, n))
get0(k, p) = get(ψ, (k, p), zeros(T, n, n))
for k = 1:4
    o = SOps(g, k * (k + 4)); H = iseven(k) ? sharmonics(g, k) : Matrix{T}[]
    wVk = lad[k][1] .* wv0 .- lad[k][2] .* wv1
    for p = (k÷2):-1:0
        wh(f) = g.w .* (2(k + 2) * (p + 1) * f + (p + 1) * (p + 2) * get0(k, p + 2)) - 2wVk .* get0(k - 1, p)
        if iseven(k) && p < k ÷ 2 && haskey(ψ, (k, p + 1))
            G = [sip(g, a, b) for a in H, b in H]
            s = [sum(g.W ./ g.w .* h .* wh(get0(k, p + 1))) for h in H]
            ψ[(k, p + 1)] += sum(c * h for (c, h) in zip(-(G \ s) ./ (2(k + 2) * (p + 1)), H))
        end
        ψ[(k, p)] = ssolve(g, o, wh(get0(k, p + 1)), H)
    end
end
xn = [cos(P * (j + T(1) / 2) / n) for j = 0:n-1]; bw = [(-1)^j * sin(P * (j + T(1) / 2) / n) for j = 0:n-1]
ax = bw ./ (T(side) .- xn); ax ./= sum(ax)
out = Dict()
for k = 1:4
    out[(0, (k, 0), 0)] = chebcoef(vec(ax' * ψ[(k, 0)]))
end
serialize(get(ENV, "OUT", "lad_edge_$(n).jls"), out)
println("saved; tails: ", [Float64(abs(out[(0, (k, 0), 0)][end])) for k = 1:4])
