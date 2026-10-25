# Null vector of the guessed ODE sum_{i<=r} sum_{j<=d} c_ij t^j B^(i)(t) = 0 (uncentered monomials), normalized.
include(joinpath(@__DIR__, "guess_ode.jl"))
L = deserialize(get(ENV, "IN", "chain_sectors.jls"))
key = eval(Meta.parse(ENV["KEY"])); r = parse(Int, ENV["R"]); d = parse(Int, ENV["D"])
c = L[key]
tlo = T(parse(Float64, get(ENV, "TLO", "0.45"))); thi = T(parse(Float64, get(ENV, "THI", "0.97")))
np = 4 * (r + 1) * (d + 1)
rows = Vector{Vector{T}}()
for i = 1:np
    t0 = tlo + (thi - tlo) * (i - 1) / (np - 1)
    s = series_at(c, t0, r)
    push!(rows, [t0^j * s[k+1] * factorial(big(k)) for k = 0:r for j = 0:d])
end
A = permutedims(hcat(rows...))
sc = [norm(A[:, k]) for k = 1:size(A, 2)]
A ./= sc'
F = qr(A, ColumnNorm())
R = F.R; m = size(R, 2)
println("diag tail: ", Float64.(abs.(diag(R))[end-2:end] ./ abs(R[1, 1])))
x = zeros(T, m); x[m] = 1
x[1:m-1] = -(UpperTriangular(R[1:m-1, 1:m-1]) \ R[1:m-1, m])
v = zeros(T, m); v[F.p] = x
v ./= sc
# normalize: leading coefficient of t^d B^(r) ... find largest abs and normalize so that coefficient (r, top nonzero) = 1
idx = findlast(x -> abs(x) > 1e-30 * maximum(abs.(v)), v)
v ./= v[idx]
for k = 0:r
    println("p_$k: ", join([string(Float64(v[k*(d+1)+j+1])) for j = 0:d], "  "))
end
open(get(ENV, "OUT", "nullvec.txt"), "w") do io
    for k = 0:r, j = 0:d
        println(io, k, " ", j, " ", string(v[k*(d+1)+j+1]))
    end
end
