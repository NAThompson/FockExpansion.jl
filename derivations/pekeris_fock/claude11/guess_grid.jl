# Coarse (r,d) grid ODE search for a far-edge/arc function (larger sizes than guess_ode.jl's exhaustive loop).
include(joinpath(@__DIR__, "guess_ode.jl"))
L = deserialize(get(ENV, "IN", "chain_sectors.jls"))
key = eval(Meta.parse(ENV["KEY"]))
c = L[key]
tlo = T(parse(Float64, get(ENV, "TLO", "1.15"))); thi = T(parse(Float64, get(ENV, "THI", "2.3")))
for r in eval(Meta.parse(get(ENV, "RS", "2:10"))), d in eval(Meta.parse(get(ENV, "DS", "4:4:40")))
    tot = (r + 1) * (d + 1)
    tot > parse(Int, get(ENV, "MAXUNK", "300")) && continue
    np = tot + 40
    rows = Vector{Vector{T}}()
    for i = 1:np
        t0 = tlo + (thi - tlo) * (i - 1) / (np - 1)
        s = series_at(c, t0, r)
        push!(rows, [((t0 - (tlo + thi) / 2) / ((thi - tlo) / 2))^j * s[k+1] * factorial(big(k)) for k = 0:r for j = 0:d])
    end
    A = permutedims(hcat(rows...))
    A ./= [norm(A[:, k]) for k = 1:size(A, 2)]'
    F = qr(A, ColumnNorm())
    dg = abs.(diag(F.R)) ./ abs(F.R[1, 1])
    @printf("r=%2d d=%2d unk=%3d  last %.1e  prev %.1e\n", r, d, tot, Float64(dg[end]), Float64(dg[end-1]))
    flush(stdout)
end
