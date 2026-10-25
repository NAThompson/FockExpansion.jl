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
function genparts(α0, θ)
    g = FE.geometry(α0, θ); π_ = T(π)
    α = FE.fold_half(g.α, π_)
    FE.primal(g.α) <= π_ / 2 || (g = merge(g, (; α, r1 = g.r2, r2 = g.r1)))
    r = sqrt(T(2)); c = cos(α); y = -g.ξ / r; h1 = g.r1 / g.r2; h2 = inv(h1); opc, omc = 2g.r1^2, 2g.r2^2
    f1(C0, C2) = FE.endpoint_quadgk((t, x) -> FE.moment_integrand(t, π_ / 2 - atan(t), atan(t / h1), h1, c, y, C0, C2, FE.quartic(c, h1, x * (1 + t), false, opc, omc), opc), -y, zero(T); rtol = rt)[1]
    f2(C0, C2) = FE.endpoint_quadgk((t, x) -> h2 * FE.moment_integrand(h2 * t, π_ / 2 - atan(t / h1), atan(t), h2, -c, y, C0, C2, FE.quartic(-c, h2, x * (1 + t), true, omc, opc), omc), -y, zero(T); rtol = rt)[1]
    z = zero(T); out = T[]
    for (f, b) in ((f1, c), (f2, -c))
        for m in (one(T), b, b^2); push!(out, f(m, z)); end
        for m in (one(T), b, y^2); push!(out, f(z, (b - 1) * m)); end
    end
    out
end
G12 = zeros(T, length(pts), 12)
for (i, p) in enumerate(pts)
    G12[i, :] = genparts(p...)
    i % 200 == 0 && (println("pt ", i); flush(stdout))
end
writedlm("genparts.txt", G12)
