include(joinpath(@__DIR__, "logcoef.jl"))
setprecision(BigFloat, 256); T = BigFloat
fs = Dict("A" => "lad_fixed_0-1101010.jls", "B" => "lad_fixed_100-11010.jls", "C" => "lad_fixed_10100-110.jls", "D" => "lad_fixed_1010100-1.jls")
B = Dict(k => deserialize(v).out[(0, 0)][(4, 0)] for (k, v) in fs)
L = sector_logcoefs(deserialize("ext_side1_n80_N60.jls"))
c1 = L[((4, 0), 1, 1)]
println("max|C1 - (A+B+C+D)| = ", Float64(maximum(abs.(c1 .+ (B["A"] + B["B"] + B["C"] + B["D"])))), "  (sign: extend uses opposite l convention? also try minus) ",
    Float64(maximum(abs.(c1 .- (B["A"] + B["B"] + B["C"] + B["D"])))))
out = Dict((0, (4, 0), k) => chebcoef(v) for (k, v) in B)
serialize("lad_B1.jls", out)
