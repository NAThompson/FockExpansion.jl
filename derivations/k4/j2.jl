include("green4.jl"); using Printf
const B=(π-2)/(3π)
ξ(a,t)=sqrt(2sin(π/4-a/2)^2+2sin(a)*sin(t/2)^2); ς(a)=cos(a/2)+sin(a/2)
V1(a)=1/sin(a/2)+1/cos(a/2)
h2d(a,t)=B*(V1(a)*(5ξ(a,t)^3/6-ξ(a,t))+ς(a)*(ξ(a,t)-1/ξ(a,t)))
for (α, θ) in ((0.7, 1.1), (1.5, 0.3), (π/2, 0.05)), (n, lv) in ((8, 8), (12, 12), (16, 16))
    t=@elapsed g=solve4(h2d, α, θ; n, levels=lv)
    @printf("(%.2f,%.2f) n=%2d levels=%2d  ψ41^(2) = % .15f  (%.2f s)\n", α, θ, n, lv, g, t)
end
