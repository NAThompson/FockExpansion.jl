# Quick Float64 test of H1: R^4 psi41^(2) in span of polynomials x order-2 atoms (+P twists), on the R = 1 slice.
using FockExpansion, LinearAlgebra, Random
using ClausenFunctions: cl2
const FE = FockExpansion

function atoms(r1, r2, r12)
    R2 = r1^2 + r2^2
    P = sqrt(2R2 - r12^2)
    d, s = r1 - r2, r1 + r2
    L5 = r12 * P + s * d
    α = 2atan(r2 / r1)
    β = asin((R2 - r12^2) / R2)
    DL = cl2(α - β) + cl2(π + β - α) - cl2(π - α - β) - cl2(α + β)
    lg = β * log((2r1 * r2)^2 / L5^2) + (β + α) * log((s^2 - r12^2) / (2r1 * r2)) + (β - α) * log((r12^2 - d^2) / (2r1 * r2))
    Ω = lg - DL
    l1, l2, l3 = log(sqrt(R2) / s), log((s + r12) / sqrt(R2)), log((d + r12)^2 / L5)
    base = [1.0, α, β, l1, l2, l3, Ω]
    vcat(base, P .* base)
end

function monos(r1, r2, r12, n)
    v = Float64[]
    for j = 0:1, i = 0:n, k = 0:n
        i + j + k <= n && push!(v, r1^i * r2^j * r12^k)
    end
    v
end

Random.seed!(1)
npts = parse(Int, get(ARGS, 1, "2500"))
pts = [(π * rand(), π * rand()) for _ = 1:npts]
vals = [FE.psi41_z2(a, t) for (a, t) in pts]
geo = [(cos(a / 2), sin(a / 2), FE.xi_stable(a, t)) for (a, t) in pts]   # R = 1; r12 = xi
for sel in (:base, :all), n in (4, 6, 8)
    rows = map(geo) do (r1, r2, r12)
        A = atoms(r1, r2, r12)
        A = sel == :base ? A[1:7] : A
        m = monos(r1, r2, r12, n)
        vcat([a .* m for a in A]...)
    end
    M = permutedims(hcat(rows...))
    F = qr(M, ColumnNorm())
    c = F \ vals
    res = norm(M * c - vals) / norm(vals)
    println("sel=$sel n=$n unknowns=$(size(M,2)) rel residual = $res")
end
