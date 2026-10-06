# Compare the stage-1 runs on different grids: relative differences of ψ_kp at fixed angles.
using DoubleFloats, Serialization, Printf
K=parse(Int, ARGS[1]); ns=parse.(Int, ARGS[2:end])
runs=Dict(n=>Dict(r[1]=>r for r in deserialize("growth_n$(n)_K$(K).jls")) for n in ns)
nref=ns[end]
kmax=minimum(maximum(keys(runs[n])) for n in ns)
for k in 2:2:kmax
    ref=runs[nref][k][3]
    sc=maximum(maximum(abs, v[2]) for v in values(ref))
    @printf("k=%3d  max|ψ_k| %.2e", k, Float64(sc))
    for n in ns[1:end-1]
        d=maximum(maximum(abs, runs[n][k][3][p][2]-ref[p][2]) for p in keys(ref))
        @printf("   n=%d vs %d: %.1e", n, nref, Float64(d/sc))
    end
    println()
end
