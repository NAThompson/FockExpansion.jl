# TASK30: independent unrolled three-step recurrence; no fit_phi3*.jl imports.
using LinearAlgebra, DelimitedFiles
include(joinpath(@__DIR__,"../../claude9/fock_spectral_copy.jl"))
BLAS.set_num_threads(1)
setprecision(BigFloat,256)
const T=BigFloat
const ROOT=normpath(joinpath(@__DIR__,"../.."))
function sample(g,Coeff,a,t)
    # Direct Chebyshev expansion, rather than Claude's barycentric evaluator.
    n=g.n; p=T(pi); beta=acos(sin(a)*cos(t)); xx=2(a-beta)/p; yy=2(a+beta-p)/p
    vals(z)=[cos(k*acos(clamp(z,-one(T),one(T)))) for k=0:n-1]
    dot(vals(xx),Coeff*vals(yy))
end
function evaluate(n,pts)
    C=[2/T(n)*cos(T(pi)*k*(j+T(1)/2)/n)/(k==0 ? 2 : 1) for k=0:n-1,j=0:n-1]
    g=SGrid(T,n); ops=[SOps(g,k*(k+4)) for k=1:3]; H=sharmonics(g,2)
    gram=[sip(g,h,j) for h in H,j in H]
    v0=2sqrt(T(2)).*sin.(g.α).*cos.(g.β./2)
    v1=4sin.(g.β).*(cos.(g.α./2).+sin.(g.α./2))
    invk(k,r)=ssolve(g,ops[k],r,k==2 ? H : Matrix{T}[])
    function path(vs)
        f1=invk(1,-2vs[1])
        proj=gram\[sum(g.W./g.w.*h.*(vs[2].*f1)) for h in H]
        f21=sum(proj[j]/4*H[j] for j in eachindex(H))
        f20=invk(2,8g.w.*f21-2vs[2].*f1)
        f31=invk(3,-2vs[3].*f21)
        f30=invk(3,10g.w.*f31-2vs[3].*f20)
        coeff=C*f30*C'
        [sample(g,coeff,p...) for p in pts]
    end
    hcat(path((v1,v0,v0)),path((v0,v1,v0)),path((v0,v0,v1)))
end

g=SGrid(T,24); a=big"0.3724764877645945"; t=big"2.29768816453775"
C=[2/T(g.n)*cos(T(pi)*k*(j+T(1)/2)/g.n)/(k==0 ? 2 : 1) for k=0:g.n-1,j=0:g.n-1]
f=g.α.^2 .+ g.β
println("poly expected=",a*a+acos(sin(a)*cos(t))," sampled=",sample(g,C*f*C',a,t))
Y=evaluate(24,[[a,t]])
println("paths=",Y)
# Coefficient from independent full-Z chain, degree <=3.
v=[]
for z in (-2,-1,1,2)
 F=sfock(g;Z=T(z),E=0,kmax=3)[(3,0)]
 push!(v,sample(g,C*F*C',a,t))
end
println("Z coefficient=",(8*(v[3]-v[2])-(v[4]-v[1]))/12," sum=",sum(Y))
