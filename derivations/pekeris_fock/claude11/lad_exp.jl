include(joinpath(@__DIR__, "logcoef.jl"))
setprecision(BigFloat, 256); T = BigFloat
fs = Dict("A" => "lad_fixed_0-1101010.jls", "B" => "lad_fixed_100-11010.jls")
out = Dict()
for (k, v) in fs, kp in ((2, 0), (3, 0))
    out[(0, kp, k)] = chebcoef(deserialize(v).out[(0, 0)][kp])
end
serialize("lad_B1_low.jls", out)
