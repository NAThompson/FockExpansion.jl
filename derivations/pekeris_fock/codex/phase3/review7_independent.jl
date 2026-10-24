# Independent exported-state integrals; includes no Claude implementation.
using LinearAlgebra, Printf, ForwardDiff
import ClausenFunctions
setprecision(BigFloat,256)
function ClausenFunctions.cl2(d::ForwardDiff.Dual{T}) where T
 t=ForwardDiff.value(d)
 ForwardDiff.Dual{T}(ClausenFunctions.cl2(t),-log(abs(2sin(t/2)))*ForwardDiff.partials(d))
end
function atom(g,x,y)
 g==:one && return one(x)
 q=(1+x*x)/2; al=big(pi)/2-2atan(x); v=1-y*y/q; P=sqrt(2q-y*y);be=asin(v)
 g==:alp && return -2atan(x)
 g==:Lt && return log(1+y)-log(q)/big(pi)
 g==:V && return P/y*(be+big(pi)/2)
 cl=ClausenFunctions.cl2; L5=y*P+x; t=(1-x*x)/2
 dl=cl(al-be)+cl(big(pi)+be-al)-cl(big(pi)-al-be)-cl(al+be)
 om=be*log(t*t/(L5*L5))+(be+al)*log((1-y*y)/t)+(be-al)*log((y*y-x*x)/t)-dl
 g==:U && return om+big(pi)*log((x+y)^2/L5)
 error(g)
end
function readstate(name)
 ls=readlines(joinpath(@__DIR__,"../../claude3/states",name));E=parse(BigFloat,split(ls[1],"E_variational = ")[2])
 ts=[(parse.(Int,split(l)[1:4])...,Symbol(split(l)[5]),parse(BigFloat,split(l)[6])) for l in ls[2:end]]
 E,ts
end
function gauss(n)
 x,v=eigen(SymTridiagonal(zeros(n),[k/sqrt(4k*k-1) for k=1:n-1])); (BigFloat.((x.+1)./2),BigFloat.(v[1,:].^2))
end
const ee=big"1.704";const gamma=BigFloat(Base.MathConstants.eulergamma)
function moment(m,p,rate)
 z=sum((inv(BigFloat(j)) for j=1:m);init=big(0))-gamma-log(rate)
 a=big(pi)^2/6-sum((inv(BigFloat(j)^2) for j=1:m);init=big(0))
 z3=5sum((-1)^(j+1)/(BigFloat(j)^3*binomial(big(2j),big(j))) for j=1:150)/2
 b=-2(z3-sum((inv(BigFloat(j)^3) for j=1:m);init=big(0)))
 c=6(big(pi)^4/90-sum((inv(BigFloat(j)^4) for j=1:m);init=big(0)))
 factorial(big(m))/rate^(m+1)*(p==0 ? one(z) : p==1 ? z : p==2 ? z*z+a : p==3 ? z^3+3z*a+b : z^4+6z*z*a+3a*a+4z*b+c)
end
function a21(nq)
 E,ts=readstate("min5_N8.txt");z,w=gauss(nq);D=big(0);nc=big(0)
 rad=Dict((n,p)=>moment(n-1,p,ee) for n=1:12 for p=0:2)
 for i=1:nq,j=1:nq
  y=z[i];x=y*z[j];q=(1+x*x)/2;v=1-y*y/q;vs=-8/(1-x*x)+1/y
  wt=64/big(pi)*y*y*(1-x*x)/(1+x*x)^3*v*w[i]*w[j];nc+=wt*v
  av=Dict(g=>atom(g,x,y) for g in unique(t[5] for t in ts));dp=big(0)
  for (n,p,a,b,g,c) in ts
   T=c*x^a*y^b*av[g]
   dp+=T*(n==0 ? vs*ee*(log(ee)+gamma-1-log(q)/2)+E*(gamma+log(ee)-log(q)/2) : n==1 ? vs*(-gamma-log(ee)+log(q)/2)-E/ee : vs*rad[(n-1,p)]-E*rad[(n,p)])
  end
  D+=wt*dp
 end
 @printf("independent a21 min5 N8 nq=%d %.15f norm=%.15f\n",nq,Float64(-D/4+2*(big(pi)-2)/(24big(pi))),Float64(nc));flush(stdout)
end
# State Rayleigh quotient assembled directly at angular nodes, with exact radial integration.
function energy(hinv)
 N=parse(Int,get(ENV,"CHECK_N","4"))
 E,ts=readstate("min5_N$(N).txt");ng=unique(t[5] for t in ts)
 # tanh-sinh tensor rule; no moment tables or basis matrices.
 zz=BigFloat[];ww=BigFloat[]
 for k=-4hinv:4hinv
  t=BigFloat(k)/hinv;u=big(pi)/2*sinh(t);z=(1+tanh(u))/2
  if big"1e-25"<z<1-big"1e-25"
   push!(zz,z);push!(ww,big(pi)/4*cosh(t)/cosh(u)^2/hinv)
  end
 end
 rm=Dict((m,p)=>moment(m,p,2ee) for m=0:2N+8 for p=0:4)
 H=big(0);S=big(0)
 for i in eachindex(zz),j in eachindex(zz)
  y=zz[i];x=y*zz[j];w=ww[i]*ww[j]*y*2; meas=(1-x*x)*y/8
  av=Dict(g=>atom(g,x,y) for g in ng)
  ag=Dict(g=>ForwardDiff.gradient(z->atom(g,z[1],z[2]),[x,y]) for g in ng)
  vals=Dict{Tuple{Int,Int},Vector{BigFloat}}()
  for (n,p,a,b,g,c) in ts
   v=c*x^a*y^b*av[g]
   dx=c*((a==0 ? big(0) : a*x^(a-1)*y^b*av[g])+x^a*y^b*ag[g][1])
   dy=c*((b==0 ? big(0) : b*x^a*y^(b-1)*av[g])+x^a*y^b*ag[g][2])
   get!(vals,(n,p),zeros(BigFloat,3)) .+= [v,dx,dy]
  end
  deriv=Dict{Tuple{Int,Int},Vector{BigFloat}}()
  for ((n,p),(v,dx,dy)) in vals
   get!(deriv,(n-1,p),zeros(BigFloat,3)) .+= [n*v+(1-x)*dx-y*dy,n*v-(1+x)*dx-y*dy,dy]
   get!(deriv,(n,p),zeros(BigFloat,3)) .+= [-ee*v,-ee*v,0]
   p>0 && (get!(deriv,(n-1,p-1),zeros(BigFloat,3)) .+= [p*v,p*v,0])
  end
  for ((n,p),a) in vals,((m,r),b) in vals
   S+=w*meas*a[1]*b[1]*rm[(n+m+5,p+r)]
   H+=w*meas*(-8/(1-x*x)+1/y)*a[1]*b[1]*rm[(n+m+4,p+r)]
  end
  A=2*(x+y*y)/((1+x)*y);B=2*(y*y-x)/((1-x)*y)
  for ((n,p),a) in deriv,((m,r),b) in deriv
   kin=a[1]*b[1]+a[2]*b[2]+2a[3]*b[3]+A*(a[1]*b[3]+a[3]*b[1])/2+B*(a[2]*b[3]+a[3]*b[2])/2
   H+=w*meas*kin*rm[(n+m+5,p+r)]/2
  end
 end
 ref=big"-2.903724377034119598311159245194404446696925310"
 @printf("independent state energy min5 N%d hinv=%d E-E0=%.15e E-Eexport=%.4e\n",N,hinv,Float64(H/S-ref),Float64(H/S-E));flush(stdout)
end
if abspath(PROGRAM_FILE)==@__FILE__
if !haskey(ENV,"CHECK_N")
 for n in (24,48,96);a21(n);end
end
for h in (haskey(ENV,"CHECK_N") ? (12,16) : (8,12));energy(h);end

end
