# Independent check: apply Λ² to ψ₄₁⁽²⁾ by finite differences and compare with
# h₄₁⁽²⁾ = 24ψ₄₂⁽²⁾ - 2V₀ψ₃₁⁽²⁾ + 2V₁ψ₃₁⁽¹⁾ (Z² part of 24ψ₄₂ - 2Vψ₃₁ + 2Eψ₂₁).
using FockExpansion, Printf
const F=FockExpansion
B=(π-2)/(3π)
function h2(a, t)
    ξ=sqrt(1-sin(a)*cos(t)); ς=cos(a/2)+sin(a/2); V1=1/sin(a/2)+1/cos(a/2)
    p31_1=B*ξ*(5ξ^2-6)/12; p31_2=B*ς*sin(a)*cos(t)/2
    24*(π-2)*(5π-14)/(180π^2)*(1-2sin(a)^2*sin(t)^2)-2p31_2/ξ+2V1*p31_1
end
ψ(a,t)=F.psi41_z2(a, t; n=10, levels=10)
for (a, t) in ((0.7, 1.1), (1.2, 2.0), (0.4, 0.5))
    hh=2e-3
    f0=ψ(a,t); fap=ψ(a+hh,t); fam=ψ(a-hh,t); ftp=ψ(a,t+hh); ftm=ψ(a,t-hh)
    faa=(fap-2*f0+fam)/hh^2; fa=(fap-fam)/(2hh); ftt=(ftp-2*f0+ftm)/hh^2; ft=(ftp-ftm)/(2hh)
    Λ2=-4*(faa+2cot(a)*fa+(ftt+cot(t)*ft)/sin(a)^2)
    @printf("(%.1f,%.1f): (Λ²-32)ψ = % .9f   h = % .9f   diff %.1e\n", a, t, Λ2-32*f0, h2(a,t), Λ2-32*f0-h2(a,t))
end
