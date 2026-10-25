# High-precision ODE search for the ladder-A far-edge variation (variable u or t), from chain.jl output files.
include(joinpath(@__DIR__, "guess_ode.jl"))   # series_at (VAR, OFF), T = BigFloat at BITS
include(joinpath(@__DIR__, "logcoef.jl"))
R1 = deserialize(ENV["F1"]); R2 = deserialize(ENV["F2"])
kp = eval(Meta.parse(get(ENV, "KP", "(4, 0)")))
c1 = chebcoef(R1.out[(0, 0)][kp]); c2 = chebcoef(R2.out[(0, 0)][kp])
err = maximum(abs(chebeval(c1, T(x)) - chebeval(c2, T(x))) for x in (-0.95, -0.6, -0.2, 0.1, 0.5, 0.9))
println("accuracy estimate |n1-n2| = ", Float64(err), "  max|B| = ", Float64(maximum(abs.(R1.out[(0, 0)][kp]))))
tlo = T(parse(Float64, get(ENV, "TLO", "-0.95"))); thi = T(parse(Float64, get(ENV, "THI", "-0.05")))
for (r, d) in eval(Meta.parse(ENV["RD"]))
    tot = (r + 1) * (d + 1); np = tot + 60
    rows = Vector{Vector{T}}()
    for i = 1:np
        t0 = tlo + (thi - tlo) * (i - 1) / (np - 1)
        s = series_at(c1, t0, r)
        push!(rows, [((t0 - (tlo + thi) / 2) / ((thi - tlo) / 2))^j * s[k+1] * factorial(big(k)) for k = 0:r for j = 0:d])
    end
    A = permutedims(hcat(rows...))
    A ./= [norm(A[:, k]) for k = 1:size(A, 2)]'
    F = qr(A, ColumnNorm())
    dg = abs.(diag(F.R)) ./ abs(F.R[1, 1])
    @printf("r=%2d d=%3d unk=%3d  last %.1e  prev %.1e  ratio %.1e\n", r, d, tot, Float64(dg[end]), Float64(dg[end-1]), Float64(dg[end] / dg[end-1]))
    flush(stdout)
end
