include("green4.jl"); using Printf
const B=(π-2)/(3π)
ξ(a,t)=sqrt(1-sin(a)*cos(t)); ς(a)=cos(a/2)+sin(a/2)
E=-2.9
h1(a,t)=(π-2)/(18π)*(6(1-2E)+(12*E-5)*ξ(a,t)^2)
p1(a,t)=(π-2)/(2880π)*(3(32*E-15)-8(12*E-5)*ξ(a,t)^2)
h3(a,t)=B*(2+tan(a/2)+cot(a/2))*sin(a)*cos(t)
p3(a,t)=-(π-2)/(120π)*(4+5sin(a))*sin(a)*cos(t)
for (α, θ) in ((0.7, 1.1), (1.5, 0.3), (2.6, 2.5))
    for (name, h, p) in (("j=1", h1, p1), ("j=3", h3, p3))
        t=@elapsed g=solve4(h, α, θ)
        @printf("%s (%.1f,%.1f): green % .12f exact % .12f diff %.1e  (%.2f s)\n", name, α, θ, g, p(α,θ), g-p(α,θ), t)
    end
end
