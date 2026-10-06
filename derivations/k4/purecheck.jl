# The per-t solutions with their degree-2 harmonic removed give a smooth integrand in
# φ (t = sin²φ), so plain Gauss–Legendre converges without endpoint grading.
include("feynman41.jl")
using FockExpansion; const F=FockExpansion
pts=((0.7,1.1),(1.5,0.3),(0.05,1.5),(2.2,2.4))
g=[F.psi41_z2_green(a,t; n=12, levels=12) for (a,t) in pts]
for n in (16, 24, 32, 48)
    r=feynman_rule(n=n, levels=0)
    tm=@elapsed v=[psi41_z2_feynman(a,t; rule=r) for (a,t) in pts]
    @printf("n=%2d (%3d nodes, %.1f ms/pt): max |Feynman - Green| = %.1e\n", n, length(r[1]), 1000tm/length(pts), maximum(abs, v-g))
end
