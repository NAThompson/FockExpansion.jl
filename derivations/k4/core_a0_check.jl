include("core_a0.jl")
using FockExpansion; const F=FockExpansion
Z=2.0
function src(a, t)
    ξ=F.xi_stable(a, t); y=-ξ/SQ2; s=0.0
    for b in (cos(a), -cos(a))
        c0=-SQ2*Z*(9Z*b^2-18Z-8b^2+8)/576
        c2=SQ2*Z*(b-1)*(90Z*y^2-45Z-18b-32y^2+16)/288
        s+=first(F.panel_adaptive(b, y, c0, c2; rtol=1e-12))
    end
    -2s/ξ
end
for (α, θ) in ((0.7, 1.1), (1.3, 2.4))
    g=F.solve_k4(src, α, θ; n=8, levels=8)
    for (nt, nu) in ((24, 24), (48, 48))
        v=core_a0(α, θ, Z; nt, nu)
        @printf("(%.1f,%.1f) nt=nu=%d: 2D zonal % .12f  Green % .12f  diff %.2e\n", α, θ, nt, v, g, v-g)
    end
end
