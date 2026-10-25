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
ids=parse.(Int,readlines(joinpath(ROOT,"codex/phase9/x19_fit40_big_ids.txt")))
lines=readlines(joinpath(ROOT,"claude9/spec40.txt"))
pts=[parse.(T,split(lines[i]))[1:2] for i in ids]
if get(ENV,"X30_MODE","eval")=="eval"
    for n in (64,80)
        Y=evaluate(n,pts)
        writedlm(joinpath(@__DIR__,"x30_targets$(n).txt"),Y)
        println("evaluated n=",n," rss=",Sys.maxrss());flush(stdout)
    end
else
    Y=readdlm(joinpath(@__DIR__,"x30_targets80.txt"),T)
    Y64=readdlm(joinpath(@__DIR__,"x30_targets64.txt"),T)
    println("spectral64_80_max=",maximum(abs.(Y-Y64)))
    Y=hcat(Y,vec(sum(Y,dims=2)))
    A=readdlm(joinpath(ROOT,"codex/phase9/x19_fit40_big_atoms.txt"),T)[:,1:115]
    Z=readdlm(joinpath(ROOT,"claude11/z1parts.txt"),T)
    G=hcat(readdlm(joinpath(ROOT,"claude11/genparts.txt"),T),readdlm(joinpath(ROOT,"claude11/genparts34.txt"),T))
    println("direct sum control max=",maximum(abs.(Y[:,4]+vec(sum(Z,dims=2)))))
    geo=[(cos(a/2),sin(a/2),sqrt(1-sin(a)*cos(t))) for (a,t) in pts]
    mons(x,y,c,d)=[x^i*y^j*c^k for j=0:1 for i=0:d for k=0:d if i+j+k<=d]
    train=get(ENV,"X30_SPLIT","original")=="original" ? collect(1:800) : [i for i=1:1000 if i%5!=0]
    test=setdiff(1:1000,train)
    cases=get(ENV,"X30_CASES","base")
    for (label,d,den) in (cases=="base" ? [("all115_linear",1,0),("all115_simplepoles",1,0)] : [("r1r2r12",3,1),("r1r2r12",4,1)])
        nm=length(mons(geo[1]...,d)); M=zeros(T,length(pts),115nm)
        for i=1:length(pts)
            m=label=="all115_simplepoles" ? T[1,inv(geo[i][1]),inv(geo[i][2]),inv(geo[i][3])] : mons(geo[i]...,d)
            for a=1:115,j=1:nm; M[i,(a-1)*nm+j]=A[i,a]*m[j];end
        end
        if label=="panels"; M=hcat(M,Z,G);end
        # Common denominator r1*r2*r12 includes individual reciprocal classes,
        # with numerator degree bound d for every classical coefficient.
        YY=den==0 ? Y : Y.*[prod(g) for g in geo]
        for j=1:size(M,2); s=norm(M[train,j]); s>0 && (M[:,j]./=s);end
        println("start ",label," d=",d," columns=",size(M,2));flush(stdout)
        F=qr(M[train,:],ColumnNorm()); diagR=abs.(diag(F.R))
        for tol in (big"1e-40",big"1e-50",big"1e-60")
            rr=findfirst(<(diagR[1]*tol),diagR); r=isnothing(rr) ? length(diagR) : rr-1
            co=F.R[1:r,1:r]\(F.Q'*YY[train,:])[1:r,:]
            residual=M[:,F.p[1:r]]*co-YY
            println(label," d=",d," tol=",tol," rank=",r," train=",[Float64(norm(residual[train,j])/norm(YY[train,j])) for j=1:4]," held=",[Float64(norm(residual[test,j])/norm(YY[test,j])) for j=1:4]);flush(stdout)
        end
        GC.gc()
    end
    println("rss=",Sys.maxrss())
end
