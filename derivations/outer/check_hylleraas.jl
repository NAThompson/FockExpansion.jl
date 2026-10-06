# Optimally truncated outer series against the Hylleraas reference at large r₂, r₁ small.
# The overall amplitude C is fitted at one point; the shape in (r₁, θ) and R is then a test.
include("outer_series.jl")
include("../radius/hylleraas.jl")
using Serialization, LinearAlgebra
setprecision(BigFloat, 512)
h=deserialize(joinpath(tempdir(), "hylleraas_12.jls"))
E=big"-2.903724377034119598311159245194404446696905"
o=outer(E, 30)
P(l, x)=l==0 ? one(x) : l==1 ? x : ((2l-1)*x*P(l-1, x)-(l-1)*P(l-2, x))/l
# ψ_out(r₁, R, cos θ) truncated at K
function ψout(r1, R, c, K)
    s=zero(r1)
    for k in 0:K, l in 0:k
        haskey(o.g, (k, l)) || continue
        s+=R^(-k)*evalg(o.g[(k, l)], r1)*P(l, c)
    end
    exp(-o.κ*R)*R^(o.σ-1)*s
end
href(r1, R, c)=hyl_eval(h, r1, R, c)
C=href(big"0.5", big"5.0", big"0.3")/ψout(big"0.5", big"5.0", big"0.3", 3)
@printf("%-18s %-10s", "(r₁, R, cosθ)", "Hylleraas")
for K in (0, 1, 2, 3, 4, 6); @printf("  K=%-8d", K); end; println()
for (r1, R, c) in ((0.5, 4, 0.3), (1.0, 5, -0.7), (0.3, 6, 0.9), (1.5, 6, 0.0), (0.8, 8, 0.5), (2.0, 8, -0.4))
    r1, R, c=BigFloat(r1), BigFloat(R), BigFloat(c)
    ref=href(r1, R, c)
    @printf("(%.1f, %d, %+.1f)      %.3e", Float64(r1), Int(R), Float64(c), Float64(ref))
    for K in (0, 1, 2, 3, 4, 6); @printf("  %.1e", Float64(abs(C*ψout(r1, R, c, K)/ref-1))); end; println()
end
