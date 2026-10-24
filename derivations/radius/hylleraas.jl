# Hylleraas reference for the helium ground state: ψ = Σ c s^l t^{2m} u^n e^{-ζs}, with one
# or more exponents ζ (a double basis converges much faster),
# s = r₁+r₂, t = r₁-r₂, u = r₁₂. With dτ ∝ u(s²-t²) ds dt du on 0 ≤ |t| ≤ u ≤ s,
#   N = ∫u(s²-t²)ψ²,
#   T = ∫[u(s²-t²)(ψ_s²+ψ_t²+ψ_u²) + 2ψ_u(s(u²-t²)ψ_s + t(s²-u²)ψ_t)],
#   V = ∫(s²-t²-4Zsu)ψ²,   E = (T+V)/N.
# All integrals reduce to I(a,b,c) = ∫ s^a t^b u^c e^{-(ζᵢ+ζⱼ)s} (b even).
using LinearAlgebra, Printf
I3(a, b, c, κ)=isodd(b) ? zero(κ) : 2/big(b+1)/(b+c+2)*factorial(big(a+b+c+2))/κ^(a+b+c+3)
# A polynomial is a Dict (a,b,c) => coefficient, times the implicit e^{-ζs}.
const Poly=Dict{NTuple{3,Int},BigFloat}
mono(a, b, c, x=big(1))=Poly((a, b, c)=>x)
function Base.:*(p::Poly, q::Poly)
    r=Poly()
    for (k1, v1) in p, (k2, v2) in q
        k=k1 .+ k2; r[k]=get(r, k, big(0))+v1*v2
    end
    r
end
Base.:+(p::Poly, q::Poly)=mergewith(+, p, q)
scale(p::Poly, x)=Poly(k=>x*v for (k, v) in p)
integrate(p::Poly, κ)=sum((v*I3(k..., κ) for (k, v) in p if all(k .>= 0)); init=zero(κ))
# basis function and its derivatives (as polynomials times e^{-ζs})
function bfun(l, m, n, ζ)
    f=mono(l, 2m, n)
    fs=scale(mono(l, 2m, n), -ζ); l>0 && (fs=fs+mono(l-1, 2m, n, big(l)))
    ft=m>0 ? mono(l, 2m-1, n, big(2m)) : Poly()
    fu=n>0 ? mono(l, 2m, n-1, big(n)) : Poly()
    (; f, fs, ft, fu)
end
function hylleraas(Ω, ζs; Z=2)
    ζs=big.(collect(ζs))
    idx=[(l, m, n, z) for z in eachindex(ζs) for l in 0:Ω for m in 0:Ω÷2 for n in 0:Ω if l+2m+n<=Ω]
    B=[bfun(i[1:3]..., ζs[i[4]]) for i in idx]
    s2t2=mono(2, 0, 0)+scale(mono(0, 2, 0), -1)
    us2t2=mono(0, 0, 1)*s2t2
    Vw=s2t2+scale(mono(1, 0, 1), -4Z)
    A1=mono(1, 0, 2)+scale(mono(1, 2, 0), -1)                  # s(u²-t²)
    A2=mono(0, 1, 0)*(mono(2, 0, 0)+scale(mono(0, 0, 2), -1))  # t(s²-u²)
    N=length(B)
    S=zeros(BigFloat, N, N); H=zeros(BigFloat, N, N)
    Threads.@threads for i in 1:N
        for j in i:N
            a, b=B[i], B[j]
            κ=ζs[idx[i][4]]+ζs[idx[j][4]]
            S[i, j]=S[j, i]=integrate(us2t2*(a.f*b.f), κ)
            kin=us2t2*(a.fs*b.fs+a.ft*b.ft+a.fu*b.fu)+A1*(a.fu*b.fs+a.fs*b.fu)+A2*(a.fu*b.ft+a.ft*b.fu)
            H[i, j]=H[j, i]=integrate(kin+Vw*(a.f*b.f), κ)
        end
    end
    L=cholesky(Symmetric(S)).L
    M=L\H/L'; M=(M+M')/2
    λ0=eigmin(Symmetric(Float64.(M)))
    F=lu(M-(big(λ0)-big(1e-10))*I); v=ones(BigFloat, N)
    for _ in 1:30; v=F\v; v/=norm(v); end
    E=dot(v, M*v)
    c=L'\v
    (; E, c, idx, ζs)
end
# ψ at (r₁, r₂, cos θ₁₂)
function hyl_eval(h, r1, r2, cθ)
    s=r1+r2; t=r1-r2; u=sqrt(r1^2+r2^2-2r1*r2*cθ)
    sum(c*s^l*t^(2m)*u^n*exp(-h.ζs[z]*s) for (c, (l, m, n, z)) in zip(h.c, h.idx))
end
