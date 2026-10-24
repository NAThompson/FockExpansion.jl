# psi41^(2) raw Feynman integral at random interior points, 256 bits, h = 1/64 (~1e-45). Output: alpha theta value (one per line).
include("psi41z2_big.jl")
using Random
npts = parse(Int, ARGS[1]); seed = parse(Int, ARGS[2]); out = ARGS[3]
Random.seed!(seed)
setprecision(BigFloat, 256)
pts = [(BigFloat(rand()) * BigFloat(π), BigFloat(rand()) * BigFloat(π)) for _ = 1:npts]
vals = Vector{BigFloat}(undef, npts)
Threads.@threads for i = 1:npts
    setprecision(BigFloat, 256) do
        vals[i] = psi41z2_raw(pts[i]...; h = BigFloat(1) / 64)
    end
end
open(out, "w") do io
    for i = 1:npts
        println(io, string(pts[i][1]), " ", string(pts[i][2]), " ", string(vals[i]))
    end
end
