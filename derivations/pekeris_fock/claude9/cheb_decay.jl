# Chebyshev coefficient decay (in X and Y separately) of psi_kp on the (X,Y) rectangle: nearest complex singularity per direction.
include("psi40_big.jl")
setprecision(BigFloat, 256)
E0 = big"-2.903724377034119598311"
n = 80
g = SGrid(BigFloat, n)
ψ = sfock(g; Z = 2, E = E0, kmax = 4, free = Dict((2, 1) => big"0.47674788"))
P = BigFloat(π)
C = [cos(P * k * (j + big(1) / 2) / n) for k = 0:n-1, j = 0:n-1] .* (2 / big(n))   # first-kind Chebyshev transform
for key in ((2, 0), (3, 0), (4, 1), (4, 0))
    A = C * ψ[key] * C'
    mx = [maximum(abs.(A[k+1, :])) for k = 0:n-1]
    my = [maximum(abs.(A[:, k+1])) for k = 0:n-1]
    sl(v, a, b) = (log10(Float64(v[b])) - log10(Float64(v[a]))) / (b - a)
    println(key, "  X-direction log10 decay per degree: ", round(sl(mx, 10, 50), digits = 3), " (10-50) ", round(sl(mx, 40, 70), digits = 3),
        " (40-70);  Y: ", round(sl(my, 10, 50), digits = 3), " ", round(sl(my, 40, 70), digits = 3),
        "   tail |c| ", Float64(mx[n]), " ", Float64(my[n]))
end
