# Checks of the continuation and dump of the far-edge log coefficients per Z-sector.
include(joinpath(@__DIR__, "logcoef.jl"))
setprecision(BigFloat, 256)
T = BigFloat
f80 = get(ENV, "F80", joinpath(@__DIR__, "ext_side1_n80_N60.jls"))
f64 = get(ENV, "F64", joinpath(@__DIR__, "ext_side1_n64_N48.jls"))
R = deserialize(f80)
P = T(π)
g = R.g; n = g.n
# (1) continuity: element-2 values vs the physical spectral interpolant just beyond the edge (tx slightly beyond side)
xn = [cos(P * (j + T(1) / 2) / n) for j = 0:n-1]
bw = [(-1)^j * sin(P * (j + T(1) / 2) / n) for j = 0:n-1]
e2 = R.e2
for z in (2,), key in ((2, 0), (3, 0), (4, 0))
    ψ = R.results[z].ψ[key]; ψ2 = R.results[z].ψ2[key]
    for i in (e2.iA - 1, e2.iA - 3)
        tx = 2 * e2.X[i] / P
        ax = bw ./ (tx .- xn); ax ./= sum(ax)
        d = maximum(abs.(vec(ax' * ψ) .- ψ2[i, :]))
        @printf("Z=%d %s  tx=%.4f  |ext - phys interp| = %.2e\n", z, string(key), Float64(tx), Float64(d))
    end
end
# (2) sector log coefficients, convergence n=64 vs n=80 (compare Chebyshev interpolants at a few ty)
L80 = sector_logcoefs(R)
R64 = deserialize(f64)
L64 = sector_logcoefs(R64)
for k in sort(collect(keys(L80)); by = string)
    c80 = chebcoef(L80[k])
    c64 = haskey(L64, k) ? chebcoef(L64[k]) : zeros(T, 64)
    d = maximum(abs(chebeval(c80, T(x)) - chebeval(c64, T(x))) for x in (-0.9, -0.5, 0.0, 0.3, 0.8))
    @printf("%-22s max|B| %.3e  |n80-n64| %.2e  tail|c| %.1e\n", string(k), Float64(maximum(abs.(L80[k]))), Float64(d), Float64(abs(c80[end])))
end
serialize(joinpath(@__DIR__, "logcoefs_side$(R.side).jls"), Dict(k => chebcoef(v) for (k, v) in L80))
