# Independent constrained Ritz implementation. No Claude includes/source/caches.
include("independent_base.jl")
setprecision(BigFloat,parse(Int,get(ENV,"BITS","192")))
const E0=big"-2.903724377034119598311159245194404446696925310"
const kap=(big(pi)-2)/(3big(pi))
function chi4(x,y)
 q=(1+x*x)/2; sn=(1-x*x)/(1+x*x);cs=2x/(1+x*x);v=1-y*y/q
 # Equivalent independently checked order-two minimal-letter reconstruction.
 h0=-(3big(pi)+10-16ClausenFunctions.cl2(big(pi)/2))/(24big(pi))
 L=log((1+y)/sqrt(q)); T=y*y/q
 (sn-2y/q)/6-v*L/3+(cs*atom(:U,x,y)+2T*atom(:V,x,y))/(6big(pi))+h0*v
end
function jets(x,y)
 q=(1+x*x)/2;d=q-y*y;p1=y/2-2
 p2=(1-2*E0)*q/12-2y/3+4*(q/3+(1-x*x)/4)-2*(big(pi)+4)/(9big(pi))*d+2q*chi4(x,y)
 b=-2kap*d
 p1,p2,b
end
# Every entry is (radial degree, log(s) power, angular value, dx,dy).
function terms(basis,x,y)
 n,a,b=basis
 if n==-1
  j=jets(x,y);jac=ForwardDiff.jacobian(z->collect(jets(z[1],z[2])),[x,y])
  p1,p2,L=j;q=(1+x*x)/2
  v2=p2+ee*p1+ee^2/2+L*log(q)/2
  dx2=jac[2,1]+ee*jac[1,1]+jac[3,1]*log(q)/2+L*x/(2q)
  dy2=jac[2,2]+ee*jac[1,2]+jac[3,2]*log(q)/2
  return [(0,0,big(1),big(0),big(0)),(1,0,p1+ee,jac[1,1],jac[1,2]),(2,0,v2,dx2,dy2),(2,1,L,jac[3,1],jac[3,2])]
 elseif n==-2
  return [(2,0,(1+x*x)/2-y*y,x,-2y)]
 end
 [(n,0,x^a*y^b,a==0 ? big(0) : a*x^(a-1)*y^b,b==0 ? big(0) : b*x^a*y^(b-1))]
end
function deriv(ts,x,y)
 out=Tuple{Int,Int,Vector{BigFloat}}[]
 for (n,p,v,dx,dy) in ts
  push!(out,(n-1,p,[n*v+(1-x)*dx-y*dy,n*v-(1+x)*dx-y*dy,dy]))
  push!(out,(n,p,[-ee*v,-ee*v,0]))
  p>0 && push!(out,(n-1,p-1,[p*v,p*v,0]))
 end
 out
end
function nodes(hinv)
 zz=BigFloat[];ww=BigFloat[]
 for k=-4hinv:4hinv
  t=BigFloat(k)/hinv;u=big(pi)/2*sinh(t);z=(1+tanh(u))/2
  if big"1e-22"<z<1-big"1e-22"
   push!(zz,z);push!(ww,big(pi)/4*cosh(t)/cosh(u)^2/hinv)
  end
 end
 zz,ww
end
function assemble(N,h)
 bas=[(-1,0,0),(-2,0,0)]
 append!(bas,[(n,a,b) for n=3:N for a=0:2:n for b=0:n-a])
 d=length(bas);H=zeros(BigFloat,d,d);S=copy(H);zz,ww=nodes(h)
 rm=Dict((m,p)=>moment(m,p,2ee) for m=0:2N+5 for p=0:2)
 println("assemble N=$N dim=$d h=$h angular_nodes=$(length(zz)^2)");flush(stdout)
 for i in eachindex(zz),j in eachindex(zz)
  y=zz[i];x=y*zz[j];wt=2y*ww[i]*ww[j]*(1-x*x)*y/8
  ts=[terms(b,x,y) for b in bas];ds=[deriv(t,x,y) for t in ts]
  V=-8/(1-x*x)+1/y;A=2*(x+y*y)/((1+x)*y);B=2*(y*y-x)/((1-x)*y)
  for k=1:d,l=1:k
   ss=big(0);hh=big(0)
   for (n,p,v,_,_) in ts[k],(m,r,z,_,_) in ts[l]
    ss+=v*z*rm[(n+m+5,p+r)]
    hh+=V*v*z*rm[(n+m+4,p+r)]
   end
   for (n,p,a) in ds[k],(m,r,b) in ds[l]
    kin=a[1]*b[1]+a[2]*b[2]+2a[3]*b[3]+A*(a[1]*b[3]+a[3]*b[1])/2+B*(a[2]*b[3]+a[3]*b[2])/2
    # zero derivative of the seed's constant has degree -1, but products have m>=3.
    hh+=kin*rm[(n+m+5,p+r)]/2
   end
   H[k,l]+=wt*hh;S[k,l]+=wt*ss
  end
 end
 for k=1:d,l=1:k-1;H[l,k]=H[k,l];S[l,k]=S[k,l];end
 bas,H,S
end
function solve4(N,h)
 bas,H,S=assemble(N,h);d=length(bas)
 scales=sqrt.(diag(S));H=H./(scales*scales');S=S./(scales*scales')
 cholesky(Symmetric(S));K=cholesky(Symmetric(H-E0*S));v=ones(BigFloat,d)
 for it=1:70
  u=K\(S*v);u/=sqrt(dot(u,S*u))
  norm(u-v)<big"1e-45" && (v=u;break)
  v=u
 end
 E=dot(v,H*v)/dot(v,S*v);res=norm(H*v-E*S*v)/(norm(H)*norm(v));c=v./scales;c/=c[1]
 println("RESULT N=$N h=$h bits=$(precision(BigFloat)) dim=$d E=$E error=$(E-E0) a21=$(c[2]) residual=$res");flush(stdout)
 open(joinpath(@__DIR__,"poly_N$(N)_h$(h).txt"),"w") do io
  println(io,"E_variational = $E")
  for (b,v) in zip(bas,c);println(io,join(b," ")," ",v);end
 end
end
function raycheck()
 for (x,y) in ((big".2",big".7"),(big".6",big".8"))
  p1,p2,L=jets(x,y); q=(1+x*x)/2
  prev=big(0)
  for s in (big"1e-4",big"1e-5",big"1e-6")
   F=sum(s^n*log(s)^p*v for (n,p,v,_,_) in terms((-1,0,0),x,y))
   target=1+s*p1+s*s*(p2+L*log(s*sqrt(q)))
   err=abs(exp(-ee*s)*F-target)
   println("ray x=$x y=$y s=$s abs_remainder=$err remainder_over_s3log=$(err/(s^3*abs(log(s))))")
   @assert err<100s^3*abs(log(s))
   prev=err
  end
 end
end
if abspath(PROGRAM_FILE)==@__FILE__
raycheck()
for N in parse.(Int,split(get(ENV,"NS","3,4,5"),","));solve4(N,parse(Int,get(ENV,"HINV","8")));end

end
