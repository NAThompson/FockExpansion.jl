include("fitbig.jl")
pts, vals = load("spec41.txt")
pp = BigFloat(π)
f(A) = [A[1], A[2], A[5], A[7] + pp * A[6], A[10] + pp / 2 * A[8]]
names = ["1", "alpha", "L", "U", "V"]
labels = String[]
for nm in names, j = 0:1, i = 0:6, k = 0:6
    i + j + k <= 6 && push!(labels, "$nm r1^$i r2^$j r12^$k")
end
geo = [(cos(a / 2), sin(a / 2), xi(a, t)) for (a, t) in pts]
M = permutedims(hcat([vcat([x .* monos(g..., 6) for x in f(atoms(g...))]...) for g in geo]...))
sc = [norm(M[:, j]) for j in 1:size(M, 2)]
c = (qr(M ./ sc') \ vals) ./ sc
# second fit on disjoint halves to estimate coefficient accuracy
c1 = (qr(M[1:1500, :] ./ sc') \ vals[1:1500]) ./ sc
c2 = (qr(M[1501:end, :] ./ sc') \ vals[1501:end]) ./ sc
open("coef5.txt", "w") do io
    for (l, x, a, b) in zip(labels, c, c1, c2)
        println(io, l, " | ", string(x), " | ", Float64(abs(a - b)))
    end
end
println("max half-split coefficient discrepancy: ", maximum(Float64.(abs.(c1 - c2))))
