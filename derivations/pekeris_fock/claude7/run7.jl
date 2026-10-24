setprecision(BigFloat, parse(Int, get(ENV, "BITS", "384")))
include("atoms7.jl")
PSI30_DIGITS[] = 50
using Serialization
const ε = big"1.704"
Nmax = parse(Int, get(ARGS, 1, "11")); variants = Symbol.(split(get(ARGS, 2, "min5split"), ","))
nd = nodes(BigFloat(1) / 16, 72); nn = length(nd.x)
atoms = Dict{Symbol,NTuple{3,Vector{BigFloat}}}(:one => (ones(BigFloat, nn), zeros(BigFloat, nn), zeros(BigFloat, nn)))
for g in (:Lt, :alp, :V, :U)
    atoms[g] = deserialize("../claude3/atomcache_$(g)_h16_k72_b384_p50.jls")
end
cache = "piececache_h16_k72_b$(precision(BigFloat))_p50.jls"
if isfile(cache)
    merge!(atoms, deserialize(cache))
else
    vals = [zeros(BigFloat, nn) for _ = 1:5]; gx = [zeros(BigFloat, nn) for _ = 1:5]; gy = [zeros(BigFloat, nn) for _ = 1:5]
    t0 = time()
    Threads.@threads for q = 1:nn
        Jm = ForwardDiff.jacobian(z -> pieces_xy(z[1], z[2]), [nd.x[q], nd.y[q]])
        v = pieces_xy(nd.x[q], nd.y[q])
        for i = 1:5
            vals[i][q] = v[i]; gx[i][q] = Jm[i, 1]; gy[i][q] = Jm[i, 2]
        end
    end
    pc = Dict(PIECES[i] => (vals[i], gx[i], gy[i]) for i = 1:5)
    serialize(cache, pc); merge!(atoms, pc)
    println("pieces tabulated ($(round(time()-t0))s)"); flush(stdout)
end
M = moments(nd, atoms, 2Nmax + 6, 2Nmax + 6)
R = Radial(ε)
open(get(ENV, "OUT", "run7.out"), "a") do io
    for N = parse(Int, get(ENV, "NMIN", "4")):Nmax, v in variants
        B = build7(N, v; ε)
        H, S = matrices(B, M, R, ε)
        E, x = ground(H, S)
        line = @sprintf("%-10s N=%2d size=%4d  E-E0=% .4e", v, N, length(B), Float64(E - E0))
        println(line); println(io, line); flush(io); flush(stdout)
    end
end
