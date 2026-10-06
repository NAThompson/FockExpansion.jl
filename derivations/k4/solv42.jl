include("s3.jl"); using Printf
Y(l, α, θ) = l==0 ? 4cos(α)^2-1 : sin(α)^2*(3cos(θ)^2-1)/2
ξ(α,θ)=sqrt(1-sin(α)*cos(θ)); ς(α)=cos(α/2)+sin(α/2)
V(α,θ,Z)=1/ξ(α,θ)-Z*(1/sin(α/2)+1/cos(α/2))
ψ21(α,θ,Z)=-Z*(π-2)/(3π)*sin(α)*cos(θ)
ψ31(α,θ,Z)=Z*(π-2)/(36π)*(6Z*ς(α)*sin(α)*cos(θ)-ξ(α,θ)*(6-5ξ(α,θ)^2))
ψ42(α,θ,Z)=(π-2)*(5π-14)/(180π^2)*Z^2*(1-2sin(α)^2*sin(θ)^2)
for n in (16, 32), Z in (1.0, 2.0), l in (0, 2)
    E=-2.9
    lhs=24sphint((a,t)->Y(l,a,t)*ψ42(a,t,Z); n)
    rhs=sphint((a,t)->Y(l,a,t)*(2V(a,t,Z)*ψ31(a,t,Z)-2E*ψ21(a,t,Z)); n)
    @printf("n=%d Z=%g l=%d  24<Y,ψ42>=% .14f   <Y,2Vψ31-2Eψ21>=% .14f\n", n, Z, l, lhs, rhs)
end
