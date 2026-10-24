include("fitbig.jl")
pts, vals = load("spec41.txt")
pts, vals = pts[1:1600], vals[1:1600]
pp = BigFloat(π)
# atoms(): 1 α β l1 l2 l3 Ω | P Pα Pβ Pl1 Pl2 Pl3 PΩ
sets = Dict(
    "S1 {1,α,l2,l3,Ω,P,Pβ}" => A -> [A[1], A[2], A[5], A[6], A[7], A[8], A[10]],
    "S2 {1,α,l2,l3,Ω,V'}" => A -> [A[1], A[2], A[5], A[6], A[7], A[10] + pp / 2 * A[8]],
    "S3 {1,α,L,Ω+πl3,V'}" => A -> [A[1], A[2], A[5], A[7] + pp * A[6], A[10] + pp / 2 * A[8]],
    "S3b {1,α,L,Ω-πl3,V'}" => A -> [A[1], A[2], A[5], A[7] - pp * A[6], A[10] + pp / 2 * A[8]],
)
for (name, f) in sort(collect(sets); by = first)
    geo = [(cos(a / 2), sin(a / 2), xi(a, t)) for (a, t) in pts]
    M = permutedims(hcat([vcat([x .* monos(g..., 6) for x in f(atoms(g...))]...) for g in geo]...))
    ntr = 1400
    sc = [norm(M[1:ntr, j]) for j in 1:size(M, 2)]
    Ms = M ./ sc'
    c = qr(Ms[1:ntr, :]) \ vals[1:ntr]
    @printf("%s unknowns=%d train %.2e test %.2e\n", name, size(M, 2), Float64(norm(Ms[1:ntr, :] * c - vals[1:ntr]) / norm(vals[1:ntr])),
        Float64(norm(Ms[ntr+1:end, :] * c - vals[ntr+1:end]) / norm(vals[ntr+1:end])))
end
