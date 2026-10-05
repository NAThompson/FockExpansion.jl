include("de_solve.jl")
Z, E, a21 = 2.0, -2.9, 0.4
function src(a, t)
    ξ=F.xi_stable(a, t); s, c=sin(a), cos(t)
    V=1/ξ-Z*(1/sin(a/2)+1/cos(a/2))
    z1=(π-2)/(2880π)*(3*(32*E-15)-8*(12*E-5)*ξ^2); z3=-(π-2)/(120π)*(4+5s)*s*c
    12*(Z*z1+Z^3*z3)-2V*psi30(a, t; Z, E, a21, rtol=1e-13)+2E*psi20(a, t; Z, E, a21)
end
h2(a, t)=12Z^2*F.psi41_z2_source(a, t)
for (α, θ) in ((0.7, 1.1), (1.3, 2.4))
    ref=psi40(α, θ; Z, E, a21, n=12, levels=12)
    for (hh, N) in ((0.25, 12), (0.125, 24), (0.0625, 48))
        tm=@elapsed v=solve_de(src, α, θ; h2, h=hh, N=N)
        @printf("(%.1f,%.1f) h=%.4f N=%2d (%5d nodes): % .13f  ref % .13f  diff %.1e  %.2f s\n", α, θ, hh, N, 4*(2N+1)^2, v, ref, v-ref, tm)
    end
end
