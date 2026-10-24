# Z^2 coefficient of psi41 (E = a21 = 0) from the spectral chain at n = 80, 256 bits, interpolated at random points.
include("psi40_big.jl")
using Random
setprecision(BigFloat, 256)
n = 80
g = SGrid(BigFloat, n)
F = [sfock(g; Z = z, E = 0, kmax = 4)[(4, 1)] for z = 1:3]
A = [big(z)^m for z = 1:3, m = 1:3]
Ai = inv(A)
F2 = Ai[2, 1] * F[1] + Ai[2, 2] * F[2] + Ai[2, 3] * F[3]
Random.seed!(parse(Int, ARGS[2]))
open(ARGS[3], "w") do io
    for _ = 1:parse(Int, ARGS[1])
        a, t = BigFloat(rand()) * BigFloat(π), BigFloat(rand()) * BigFloat(π)
        println(io, string(a), " ", string(t), " ", string(interp(g, F2, a, t)))
    end
end
