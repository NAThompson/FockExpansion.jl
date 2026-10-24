# Reference data from the 256-bit spectral chain (fock_spectral.jl):
#  * the pure-Z harmonic coefficients of ψ₄₁ (PSI41_HARMONIC in src/fourth_order.jl), and
#  * point values of ψ₄₁ and ψ₄₀ for test/fourth_order.jl.
# Two grid sizes give the accuracy.
include("fock_spectral.jl")
setprecision(BigFloat, 256)
# Coefficients of Y₄₀ = 4cos²α-1 and Y₄₂ = sin²α P₂(cos θ) in ψ₄₁.
function harmonic41(g, ψ)
    c=hcoef(HBasis(g, 4), ψ[(4, 1)])
    c[1], c[3]
end
ns=(48, 64)
grids=[SGrid(BigFloat, n) for n in ns]
# ψ₄₁'s harmonic part is Σ_m Z^m (c0_m Y₄₀ + c2_m Y₄₂) at E = a₂₁ = 0, m = 1..3.
coef=map(grids) do g
    rows=[harmonic41(g, sfock(g; Z=big(z), E=0, kmax=4)) for z in 1:3]
    A=[big(z)^m for z in 1:3, m in 1:3]
    A\first.(rows), A\last.(rows)
end
for (name, i) in (("l0", 1), ("l2", 2))
    for m in 1:3
        a, b=coef[1][i][m], coef[2][i][m]
        @printf("%s Z^%d: %s  (n=48 vs 64: %.1e)\n", name, m, string(round(b; sigdigits=40)), Float64(abs(a-b)))
    end
end
# Point values with generic constants.
Z, E, a21, a40, a42=big"1.7", big"-2.2", big"0.31", big"0.1", big"-0.2"
pts=((big"0.7", big"1.1"), (big"1.5", big"0.3"), (big"2.3", big"2.0"))
vals=map(grids) do g
    ψ=sfock(g; Z, E, kmax=4, free=Dict((2, 1)=>a21, (4, 0)=>a40, (4, 2)=>a42))
    [(interp(g, ψ[(4, 1)], a, t), interp(g, ψ[(4, 0)], a, t)) for (a, t) in pts]
end
for (j, (a, t)) in enumerate(pts), (name, i) in (("psi41", 1), ("psi40", 2))
    v1, v2=vals[1][j][i], vals[2][j][i]
    @printf("%s(%s, %s) = %s  (n=48 vs 64: %.1e)\n", name, string(a), string(t), string(round(v2; sigdigits=25)), Float64(abs(v1-v2)))
end
