# RHS of the direct-operator check, evaluated from stencils regenerated in memory (no stored stencil files).
# Input: terms file, arrays file ("k|atom|i j q|num|den"), points file ("x y z" rationals as num/den strings), eps = 3/2.
# Output: lines "point_index|k|target_atom|G" with G = sum over blocks (sigma2P, C1, C2 with weights) of sum_e (rec A)(e) L_e(U)L_e(V)L_e(W).
include("stencil3.jl")
setprecision(BigFloat, 160)
blocks = read_terms(ARGS[1])
A = Dict{Tuple{Int,String},Dict{NTuple{3,Int},Rational{BigInt}}}()
for line in eachline(ARGS[2])
    k, at, idx, num, den = split(line, '|')
    d = get!(A, (parse(Int, k), at), Dict{NTuple{3,Int},Rational{BigInt}}())
    d[Tuple(parse.(Int, split(idx)))] = parse(BigInt, num) // parse(BigInt, den)
end
pts = [Tuple(parse(BigFloat, s) for s in split(l)) for l in eachline(ARGS[3])]
const EPSB = big(3) / 2
lag(n, t) = (n == 0 ? one(t) : (p0 = one(t); p1 = 1 - t; for k = 1:n-1; p0, p1 = p1, ((2k + 1 - t) * p1 - k * p0) / (k + 1); end; p1))
evalc(R, e) = sum(BigFloat(c) * BigFloat(e[1])^mk[1] * BigFloat(e[2])^mk[2] * BigFloat(e[3])^mk[3] * EPSB^mk[4] * BigFloat(π)^mk[5] for (mk, c) in R)
G = Dict{Tuple{Int,Int,String},BigFloat}()
for key in sort(collect(keys(blocks)))
    name, tgt, src = Meta.parse(replace(key, "'" => "\"")) |> eval
    tgt = split(tgt, ':')[2]; src = split(src, ':')[2]
    for k in (0, 1)
        kk, w8 = name == "sigma2P" ? (k, 1) : name == "C1" ? (k + 1, k + 1) : (k + 2, (k + 1) * (k + 2))
        haskey(A, (kk, src)) || continue
        arr = A[(kk, src)]
        rec = block_recurrence(blocks[key])
        for (pi_, (x, y, z)) in enumerate(pts)
            U, V, W = EPSB * (y + z - x), EPSB * (x + z - y), 2EPSB * (x + y - z)
            g = zero(BigFloat)
            Lc = Dict{NTuple{3,Int},BigFloat}()
            for (sidx, aval) in arr, (o, R) in rec
                e = (sidx[1] - o[1], sidx[2] - o[2], sidx[3] - o[3])
                minimum(e) < 0 && continue
                Le = get!(() -> lag(e[1], U) * lag(e[2], V) * lag(e[3], W), Lc, e)
                g += w8 * BigFloat(aval) * evalc(R, e) * Le
            end
            G[(pi_, k, tgt)] = get(G, (pi_, k, tgt), zero(BigFloat)) + g
        end
    end
end
open(ARGS[4], "w") do io
    for ((p, k, t), g) in G
        println(io, p, "|", k, "|", t, "|", string(g))
    end
end
println("done ", length(G))
