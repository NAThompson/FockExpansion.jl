include(joinpath(@__DIR__, "logcoef.jl"))
setprecision(BigFloat, 256); T = BigFloat
L = deserialize("lad_B1.jls"); L2 = deserialize("lad_B1_low.jls"); S = deserialize("logcoefs_side1.jls")
for (nm, c) in [("A4", L[(0, (4, 0), "A")]), ("B4", L[(0, (4, 0), "B")]), ("D4", L[(0, (4, 0), "D")]), ("A3", L2[(0, (3, 0), "A")]), ("B3", L2[(0, (3, 0), "B")]),
                ("C1", S[((4, 0), 1, 1)]), ("C2", S[((4, 0), 1, 2)]), ("C3", S[((4, 0), 1, 3)]), ("psi30Z1", S[((3, 0), 1, 1)])]
    println(nm, " e2-end(ty=1): ", string(round(chebeval(c, T(1)), sigdigits = 55)), "   -e1-end(ty=-1): ", string(round(chebeval(c, T(-1)), sigdigits = 55)))
end
