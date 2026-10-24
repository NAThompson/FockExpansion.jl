# Restrict Z-sectors of psi30/psi40 (E = a21 = 0) to the collinear edge X = -pi/2 (theta = pi, alpha in (0, pi/2)),
# parametrised by tau = tan(alpha/2) = r2/r1, and write f^(i)(tau_m)/i!, i = 0..K, at sample points tau_m.
include("psi40_big.jl")
bits = parse(Int, get(ENV, "BITS", "256")); n = parse(Int, get(ENV, "NGRID", "80")); K = parse(Int, get(ENV, "KORD", "14"))
setprecision(BigFloat, bits)
T = BigFloat
g = SGrid(T, n)
P = T(π)
# truncated power series arithmetic (coefficients in ε)
smul(a, b) = [sum(a[i+1] * b[k-i+1] for i = 0:k) for k = 0:length(a)-1]
function sinv(a)
    b = zeros(T, length(a)); b[1] = 1 / a[1]
    for k = 1:length(a)-1
        b[k+1] = -sum(a[i+1] * b[k-i+1] for i = 1:k) / a[1]
    end
    b
end
function sint(a, c0)   # integrate series: returns c0 + ∫
    vcat([c0], [a[k] / k for k = 1:length(a)-1])
end
function edge_fun(F)
    xn = [cos(P * (j + T(1) / 2) / n) for j = 0:n-1]
    bw = [(-1)^j * sin(P * (j + T(1) / 2) / n) for j = 0:n-1]
    ax = bw ./ (-1 .- xn); ax ./= sum(ax)
    h = vec(ax' * F)                         # values at Y nodes (ty = xn)
    # Chebyshev coefficients of h(ty)
    [2 / T(n) * sum(h[j+1] * cos(P * k * (j + T(1) / 2) / n) for j = 0:n-1) / (k == 0 ? 2 : 1) for k = 0:n-1]
end
function series_at(c, τ0)
    # ty(τ) = 8 atan(τ)/π - 1 ; atan series via ∫ 1/(1+τ²)
    e = zeros(T, K + 1); e[1] = τ0; e[2] = 1
    d = sinv(smul(e, e) .+ [1; zeros(T, K)])
    at = sint(d, atan(τ0))
    ty = 8 .* at ./ P; ty[1] -= 1
    # Clenshaw with series
    b1 = zeros(T, K + 1); b2 = zeros(T, K + 1)
    for k = length(c)-1:-1:1
        b1, b2 = 2 .* smul(ty, b1) .- b2 .+ [c[k+1]; zeros(T, K)], b1
    end
    smul(ty, b1) .- b2 .+ [c[1]; zeros(T, K)]
end
sectors = Dict{String,Matrix{T}}()
F3 = [sfock(g; Z = z, E = 0, kmax = 4) for z = 0:4]
A3 = inv([big(z)^m for z = 0:3, m = 0:3]); A4 = inv([big(z)^m for z = 0:4, m = 0:4])
for j = 0:3
    sectors["psi30_Z$j"] = sum(A3[j+1, i] * F3[i][(3, 0)] for i = 1:4)
end
for j = 0:4
    sectors["psi40_Z$j"] = sum(A4[j+1, i] * F3[i][(4, 0)] for i = 1:5)
end
sectors["psi20_Z1"] = sum(A3[2, i] * F3[i][(2, 0)] for i = 1:4)
sectors["psi41_Z2"] = sum(A3[3, i] * F3[i][(4, 1)] for i = 1:4)
taus = [T(i) / 200 for i = 20:180]      # tau in [0.1, 0.9]
open("edge_series.txt", "w") do io
    for (name, F) in sort(collect(sectors); by = first)
        c = edge_fun(F)
        for τ in taus
            println(io, name, " ", string(τ), " ", join(string.(series_at(c, τ)), " "))
        end
    end
end
println("done; tail Chebyshev |c| of psi40_Z3 edge: ", Float64(abs(edge_fun(sectors["psi40_Z3"])[end])))
