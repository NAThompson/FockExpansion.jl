include(joinpath(@__DIR__, "logcoef.jl"))
setprecision(BigFloat, 576); T = BigFloat
for f in ("ladA_hi_160.jls", "ladA_hi_136.jls")
    R = deserialize(f)
    for kp in ((3, 0), (4, 0))
        c = chebcoef(R.out[(0, 0)][kp])
        println(f, " ", kp, " e2-end: ", string(round(chebeval(c, T(1)), sigdigits = 100)))
    end
end
