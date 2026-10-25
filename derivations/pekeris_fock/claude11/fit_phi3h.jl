# Control psi30^(1) = classical^(1) + elliptic^(1) (repo closed forms), then identify the ladder pieces phi3^A, phi3^B.
using FockExpansion, LinearAlgebra, DelimitedFiles
const FE = FockExpansion
include(joinpath(@__DIR__, "fock_spectral_import.jl"))
setprecision(BigFloat, 256); T = BigFloat
X19 = joinpath(@__DIR__, "..", "codex", "phase9")
ids = parse.(Int, readlines(joinpath(X19, "x19_fit40_big_ids.txt")))[1:parse(Int, get(ENV, "NPTS", "1000"))]
lines = readlines(joinpath(@__DIR__, "..", "claude9", "spec40.txt"))
pts = [parse.(T, split(lines[i]))[1:2] for i in ids]
A = readdlm(joinpath(X19, "x19_fit40_big_atoms.txt"), BigFloat)[1:length(ids), :]
geometryx(a, t) = (cos(a / 2), sin(a / 2), sqrt(2sin(T(π) / 4 - a / 2)^2 + 2sin(a) * sin(t / 2)^2))
geo = [geometryx(p...) for p in pts]
rt = T(10)^-55
function genparts34(α0, θ)
    g = FE.geometry(α0, θ); π_ = T(π)
    α = FE.fold_half(g.α, π_)
    FE.primal(g.α) <= π_ / 2 || (g = merge(g, (; α, r1 = g.r2, r2 = g.r1)))
    b = -g.v; opb, omb = g.ξ^2, g.opv; h = sqrt(opb / omb); uh = FE.primal(h) <= 1; up = uh ? h : one(h)
    f(yy, C0, C2) = FE.endpoint_quadgk((t, x) -> (q = up * t; up * FE.moment_integrand(q, π_ / 2 - atan(q), atan(q / h), h, b, yy, C0, C2, FE.quartic(b, h, x * (1 + t), uh, opb, omb), opb)), min(g.r1, g.r2), zero(T); rtol = rt)[1]
    z = zero(T); out = T[]
    for yy in (-g.r1, -g.r2)
        for m in (opb * omb, one(T), b); push!(out, f(yy, m, z)); end
        for m in (one(T), yy^2, b); push!(out, f(yy, z, (b - 1) * m)); end
    end
    out
end
G12 = zeros(T, length(pts), 12)
for (i, p) in enumerate(pts)
    G12[i, :] = genparts34(p...)
    i % 200 == 0 && (println("pt ", i); flush(stdout))
end
writedlm("genparts34.txt", G12)
