# PSI41_Z2_SHIFT to high precision. The harmonic correction (2/π)M(ρ)U₂(p̂·x) removed
# from each per-t zonal solution integrates to b₀Y₄₀ + b₂Y₄₂ (exactly a k = 4 harmonic),
# so evaluating its t-integral at two points determines b. M(ρ) is closed form; the
# integrand grows like log t at t → 0, handled by adaptive BigFloat quadrature.
include("feynman41.jl")
using FockExpansion; const F=FockExpansion
Ucorr(kind, p1, p2, ω, x)=begin
    ρ=sqrt(p1^2+p2^2)
    c=(p1*x[1]+p2*x[2])/ρ
    (2/oftype(ρ, π))*Mzonal(kind, ρ, ω)*(4c^2-1)
end
function corr_integrand(t, omt, x, σ)
    Flo(v, kind)=Ucorr(kind, v[1], 1-v[2], 2v[2]-v[2]^2-v[1]^2, x)
    Fhi(v, kind)=Ucorr(kind, σ*(1-v[1]), v[2], 2v[1]-v[1]^2-v[2]^2, x)
    F_, v0, J=t<=1/2 ? (Flo, [σ*t, t], [1 0; 0 -1]) : (Fhi, [omt, omt], [-σ 0; 0 1])
    J=convert(Matrix{typeof(t)}, J)
    d1=J*ForwardDiff.gradient(v->F_(v, :log), v0)
    H=J*ForwardDiff.hessian(v->F_(v, :xlog), v0)*J
    -F_(v0, :inv)/3+(7/3)*d1[2]+(5/3)*H[2, 2]+σ*H[1, 2]
end
setprecision(BigFloat, 320)
B=(big(π)-2)/(3big(π))
function corr_total(α, θ)
    x=BigFloat[cos(α), sin(α)*cos(θ), sin(α)*sin(θ), 0]
    f(φ)=(t=sin(φ)^2; omt=cos(φ)^2; 2*(corr_integrand(t, omt, x, 1)+corr_integrand(t, omt, x, -1)))
    B/(sqrt(big(2))*big(π))*quadgk(f, big(0), big(π)/4, big(π)/2; rtol=big(10)^-40)[1]
end
Y40(α,θ)=4cos(α)^2-1; Y42(α,θ)=sin(α)^2*(3cos(θ)^2-1)/2
pts=((big"0.7", big"1.1"), (big"1.3", big"0.4"), (big"0.4", big"2.2"))
d=[corr_total(a, t) for (a, t) in pts]
M=[Y40(a,t) for (a,t) in pts]; M=hcat(M, [Y42(a,t) for (a,t) in pts])
b=M[1:2, :]\d[1:2]
println("b₀ = ", round(b[1], sigdigits=35))
println("b₂ = ", round(b[2], sigdigits=35))
println("check at third point: ", Float64(d[3]-M[3,:]'*b))
