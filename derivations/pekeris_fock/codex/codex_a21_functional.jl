# Independent radial Mellin finite parts; no Claude source is included.
using LinearAlgebra, Printf
function gauss(n)
    x,V=eigen(SymTridiagonal(zeros(n),[k/sqrt(4k*k-1) for k=1:n-1]))
    ((x .+ 1)./2, V[1,:].^2)
end
function estimate(path,n; bits=128)
 setprecision(BigFloat,bits) do
    lines=readlines(path); E=parse(BigFloat,split(lines[1],"E_variational = ")[2]); e=BigFloat("1.704")
    terms=[(parse.(Int,split(l)[1:4]),parse(BigFloat,split(l)[5])) for l in lines[2:end]]
    xq,wq=gauss(n); D=BigFloat(0); normcheck=BigFloat(0); localval=BigFloat(0)
    gamma=BigFloat(Base.MathConstants.eulergamma); c=log(2e)
    # integral s^(h-1) exp(-e*s) (log(s)+c)^p
    L(h,p)=begin
        value=BigFloat(factorial(big(h-1)))/e^h
        p==0 ? value : value*(sum(inv(BigFloat(j)) for j=1:h-1;init=BigFloat(0))-gamma-log(e)+c)
    end
    radial = [L(h,p) for h=1:13,p=0:1]
    for (iy,yf) in enumerate(xq),(it,tf) in enumerate(xq)
        y=BigFloat(yf); x=y*BigFloat(tf); hx=sqrt((1+x*x)/2)
        vv=1-2y*y/(1+x*x); Vs=-8/(1-x*x)+1/y
        weight=64/BigFloat(pi)*y*(1-x*x)/(1+x*x)^3*vv*y*BigFloat(wq[iy])*BigFloat(wq[it])
        normcheck+=weight*vv
        U=e*(y-x); V=e*(y+x); W=2e*(1-y)
        powers=[ [z^j for j=0:13] for z in (U,V,W)]
        dp=BigFloat(0); p0=BigFloat(0); p1=BigFloat(0); p2=BigFloat(0); l2=BigFloat(0)
        for (key,coef) in terms
            a,b,d,p=key; k=a+b+d
            T=coef*powers[1][a+1]*powers[2][b+1]*powers[3][d+1]
            if k==0
                @assert p==0
                dp+=T*(Vs*e*(log(e)+gamma-1-log(hx))+E*(gamma+log(e)-log(hx)))
            elseif k==1
                @assert p==0
                dp+=T*(Vs*(-gamma-log(e)+log(hx))-E/e)
            else
                dp+=T*(Vs*radial[k-1,p+1]-E*radial[k,p+1])
            end
            if p==0
                k==0 && (p0+=T); k==1 && (p1+=T); k==2 && (p2+=T)
            elseif p==1 && k==2
                l2+=T
            end
        end
        D+=weight*dp
        localval+=weight*(p2-e*p1+e^2*p0/2+log(2e/hx)*l2)/hx^2
    end
    b=-2*(BigFloat(pi)-2)/(3BigFloat(pi))
    (-D/4-b/8,localval,normcheck)
 end
end
# Manufactured f=R^2 log R exp(-R): constant Mellin source = -1.
# Gamma'(z+2)-9Gamma'(z+1)+8Gamma(z)-2Gamma(z+1).
gamma=BigFloat(Base.MathConstants.eulergamma)
fp=(1-gamma)+9gamma-8gamma-2
@assert abs(fp+1)<big"1e-60"
@assert abs(-fp/8-BigFloat(1)/8)<big"1e-60"
println("Manufactured radial check: FP=-1, a=0 PASS")
for n in (24,48,96)
 for N=6:13
    a,localv,normv=estimate(joinpath(@__DIR__,"../julia/states/fixedlog_N$(N).txt"),n)
    @printf("n=%d N=%d global=%.12f local=%.12f P[v]=%.12f\n",n,N,Float64(a),Float64(localv),Float64(normv)); flush(stdout)
 end
end
