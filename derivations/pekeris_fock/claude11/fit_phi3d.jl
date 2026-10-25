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
# panel pieces, Z^1 parts, separated: (panel1 c0, panel1 c2, panel2 c0, panel2 c2)
function z1parts(α0, θ)
    g = FE.geometry(α0, θ); π_ = T(π)
    α = FE.fold_half(g.α, π_)
    FE.primal(g.α) <= π_ / 2 || (g = merge(g, (; α, r1 = g.r2, r2 = g.r1)))
    r = sqrt(T(2)); c = cos(α); y = -g.ξ / r; h1 = g.r1 / g.r2; h2 = inv(h1); opc, omc = 2g.r1^2, 2g.r2^2
    c0 = -r / 576 * (-8c^2 + 8)
    c2a = r / 288 * (c - 1) * (-18c - 32y^2 + 16); c2b = r / 288 * (-c - 1) * (18c - 32y^2 + 16)
    out = T[]
    for (C0, C2a, C2b) in ((c0, 0, 0), (0, c2a, 0), (0, 0, 0), (0, 0, 0))
    end
    f1(C0, C2) = FE.endpoint_quadgk((t, x) -> FE.moment_integrand(t, π_ / 2 - atan(t), atan(t / h1), h1, c, y, C0, C2, FE.quartic(c, h1, x * (1 + t), false, opc, omc), opc), -y, zero(T); rtol = rt)[1]
    f2(C0, C2) = FE.endpoint_quadgk((t, x) -> h2 * FE.moment_integrand(h2 * t, π_ / 2 - atan(t / h1), atan(t), h2, -c, y, C0, C2, FE.quartic(-c, h2, x * (1 + t), true, omc, opc), omc), -y, zero(T); rtol = rt)[1]
    [f1(c0, zero(T)), f1(zero(T), c2a), f2(c0, zero(T)), f2(zero(T), c2b)], α, g
end
cl1 = T[]; P4 = zeros(T, length(pts), 4)
for (i, p) in enumerate(pts)
    v, α, g = z1parts(p...)
    P4[i, :] = v
    p0 = 2atan(g.ξ, sqrt(g.opv))
    push!(cl1, (FE.classical(α, p0, 1) - FE.classical(α, p0, -1)) / 2)
    i % 100 == 0 && (println("pt ", i); flush(stdout))
end
writedlm("z1parts.txt", hcat(cl1, P4))
