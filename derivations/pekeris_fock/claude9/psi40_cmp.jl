# BigFloat psi40 via the repo's precision-generic spectral chain (copied: fock_spectral_copy.jl), checked against repo Float64 psi40.
include("fock_spectral_copy.jl")
using FockExpansion
function interp(g::SGrid{T}, F, α, θ) where {T}
    P = T(π)
    β = acos(sin(α) * cos(θ))
    tx, ty = 2(α - β) / P, 2(α + β - P) / P
    x = [cos(P * (j + T(1) / 2) / g.n) for j = 0:g.n-1]
    bw = [(-1)^j * sin(P * (j + T(1) / 2) / g.n) for j = 0:g.n-1]
    ax = bw ./ (tx .- x)
    ay = bw ./ (ty .- x)
    (ax' * F * ay) / (sum(ax) * sum(ay))
end
if abspath(PROGRAM_FILE) == @__FILE__
    setprecision(BigFloat, 256)
    E0 = big"-2.903724377034119598311"
    a21 = big"0.47674788"     # any value; compared consistently
    pts = ((0.7, 1.1), (1.5, 0.3), (2.3, 2.0), (0.2, 2.9), (1.0, 3.1), (0.6, 3.14159), (1.5707963267948966, 1.3))
    res = Dict()
    for n in (64, 80)
        t = @elapsed begin
            g = SGrid(BigFloat, n)
            ψ = sfock(g; Z = 2, E = E0, kmax = 4, free = Dict((2, 1) => a21))
        end
        res[n] = [interp(g, ψ[(4, 0)], big(a), big(t)) for (a, t) in pts]
        println("n=$n  time $(round(t, digits=1)) s")
    end
    for (i, (a, t)) in enumerate(pts)
        r = FockExpansion.psi40(a, t; Z = 2.0, E = Float64(E0), a21 = Float64(a21))
        println((a, t), " psi40 = ", round(res[80][i], sigdigits = 45), "  |64-80| = ", Float64(abs(res[64][i] - res[80][i])),
            "  repo-Float64 diff = ", r - Float64(res[80][i]))
    end
end
