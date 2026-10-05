# Fit PSI41_Z2_SHIFT: the Y₄₀, Y₄₂ components of the Feynman sum, against the Green's
# function solution (which has none), and check the fixed rule's accuracy.
using FockExpansion, Printf
const F=FockExpansion
include(joinpath(@__DIR__, "feynman41.jl"))
pts=[(0.7,1.1),(1.5,0.3),(1.2,2.0),(0.3,0.5),(2.6,2.5),(1.0,0.05),(π/2,3.0),(0.05,1.5),(0.9,2.9),(1.3,1.6)]
g=[F.psi41_z2_green(a,t; n=12, levels=12) for (a,t) in pts]
M=hcat([F.harmonic40(a,t) for (a,t) in pts], [F.harmonic42(a,t) for (a,t) in pts])
fine=feynman_rule(n=20, levels=24)
vf=[psi41_z2_feynman(a,t; rule=fine) for (a,t) in pts]
b=M\(vf-g)
@printf("b = (%.16f, %.16f), fit residual %.1e\n", b..., maximum(abs, vf-g-M*b))
for r in (feynman_rule(n=8, levels=10), feynman_rule(n=10, levels=12), feynman_rule())
    v=[psi41_z2_feynman(a,t; rule=r) for (a,t) in pts]
    @printf("rule with %d nodes: max |Δ| vs fine rule %.1e\n", length(r[1]), maximum(abs, v-vf))
end
for (n, lv) in ((14, 16), (16, 18), (12, 18), (16, 16))
    r=feynman_rule(n=n, levels=lv)
    v=[psi41_z2_feynman(a,t; rule=r) for (a,t) in pts]
    @printf("n=%d levels=%d (%d nodes): max |Δ| vs fine rule %.1e\n", n, lv, length(r[1]), maximum(abs, v-vf))
end
