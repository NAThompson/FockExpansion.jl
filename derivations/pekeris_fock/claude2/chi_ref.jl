# reference values of the repo's chi at (r1, r2, r12) points (R = 1), both orderings of r1, r2
using FockExpansion
setprecision(BigFloat, 256)
pts = [(0.3, 0.7, 0.6), (0.8, 0.2, 0.75), (0.55, 0.45, 0.2), (0.9, 0.15, 0.95), (0.1, 0.6, 0.55)]
for (a, b, c) in pts
    R = sqrt(big(a)^2 + big(b)^2)
    r1, r2, r12 = big(a) / R, big(b) / R, big(c) / R
    α = 2atan(r2 / r1)
    μ = (r1^2 + r2^2 - r12^2) / (2r1 * r2)
    println(a, " ", b, " ", c, " ", FockExpansion.chi(α, μ))
end
