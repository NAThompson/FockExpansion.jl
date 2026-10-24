# Charge sectors C_j (Z^j coefficients, j = 0..4) of pure psi40 at E = a21 = a40 = a42 = 0, spectral chain n = 80, 256 bits.
# Output per line: alpha theta C0 C1 C2 C3 C4
include("psi40_big.jl")
using Random
setprecision(BigFloat, 256)
n = 80
g = SGrid(BigFloat, n)
F = [sfock(g; Z = z, E = 0, kmax = 4)[(4, 0)] for z = 0:4]
A = [big(z)^m for z = 0:4, m = 0:4]
Ai = inv(A)
Cs = [sum(Ai[j+1, i] * F[i] for i = 1:5) for j = 0:4]
Random.seed!(parse(Int, ARGS[2]))
open(ARGS[3], "w") do io
    for _ = 1:parse(Int, ARGS[1])
        a, t = BigFloat(rand()) * BigFloat(π), BigFloat(rand()) * BigFloat(π)
        println(io, join(string.([a, t, [interp(g, C, a, t) for C in Cs]...]), " "))
    end
end
