# Z-sector iterated variations N^s psi_kp at the far end of step s, n=80 vs n=64 convergence, Chebyshev series dump.
using Serialization, LinearAlgebra, Printf
include(joinpath(@__DIR__, "logcoef.jl"))   # chebcoef, chebeval
setprecision(BigFloat, 256)
T = BigFloat
function sectors(R)
    Zs = sort(unique(first.(collect(keys(R.out)))))
    V = inv([T(z)^m for z in Zs, m = 0:length(Zs)-1])
    d = Dict{Any,Vector{T}}()
    for s = 0:R.steps, kp in keys(R.out[(Zs[1], s)]), j = 0:length(Zs)-1
        v = sum(V[j+1, i] * R.out[(Zs[i], s)][kp] for i in eachindex(Zs))
        maximum(abs.(v)) > T(10)^-30 && (d[(s, kp, j)] = v)
    end
    d
end
R80 = deserialize(get(ENV, "F80", joinpath(@__DIR__, "chain_n80_N60.jls")))
R64 = deserialize(get(ENV, "F64", joinpath(@__DIR__, "chain_n64_N48.jls")))
S80 = sectors(R80); S64 = sectors(R64)
for k in sort(collect(keys(S80)); by = string)
    c80 = chebcoef(S80[k]); c64 = haskey(S64, k) ? chebcoef(S64[k]) : zeros(T, 64)
    dd = maximum(abs(chebeval(c80, T(x)) - chebeval(c64, T(x))) for x in (-0.9, -0.5, 0.0, 0.3, 0.8))
    @printf("step %d %-8s Z^%d  max %.3e  |n80-n64| %.2e  tail %.1e\n", k[1], string(k[2]), k[3], Float64(maximum(abs.(S80[k]))), Float64(dd), Float64(abs(c80[end])))
end
# also list sector entries present at n=64 but absent at n=80 (should be none)
for k in keys(S64)
    haskey(S80, k) || println("only in n=64: ", k, " max ", Float64(maximum(abs.(S64[k]))))
end
serialize(joinpath(@__DIR__, "chain_sectors.jls"), Dict(k => chebcoef(v) for (k, v) in S80))
