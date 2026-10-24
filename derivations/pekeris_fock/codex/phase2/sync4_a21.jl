# Independent extension of Codex's Mellin functional; no Claude source included.
using LinearAlgebra, Printf
import ClausenFunctions
function gauss(n)
 x,V=eigen(SymTridiagonal(zeros(n),[k/sqrt(4k*k-1) for k=1:n-1])); ((x.+1)./2,V[1,:].^2)
end
function chi(x,y)
 a=(1+x)/2;b=(1-x)/2;d=y;q=a*a+b*b
 al=2atan(b/a);sn=2a*b/q;cs=(a*a-b*b)/q;v=(q-d*d)/q
 gg=d*sqrt(2q-d*d)/q;sig=(a+b)/sqrt(q);xi=d/sqrt(q);be=asin(v)
 cl(t)=ClausenFunctions.cl2(t)
 dl=cl(al-be)+cl(big(pi)+be-al)-cl(big(pi)-al-be)-cl(al+be)
 logs=be*(2log(sn)-2log(gg+cs))+(be+al)*log1p(v/sn)+(be-al)*log1p(-v/sn)
 h0=-(3big(pi)+10-16BigFloat(Base.MathConstants.catalan))/(24big(pi))
 (gg+sn-2sig*xi)/6+(-2big(pi)*v*log(sig+xi)+big(pi)*cs*log((cs+sig*xi)^2/(sig^2*(gg+cs)))+2gg*be+cs*(logs-dl))/(6big(pi))+h0*v
end
function states(path)
 lines=readlines(path); E=parse(BigFloat,split(lines[1],"E_variational = ")[2])
 terms=[(parse.(Int,split(l)[1:4]),Symbol(split(l)[5]),parse(BigFloat,split(l)[6])) for l in lines[2:end]]
 E,terms
end
function run(nq,bits,Ns)
 setprecision(BigFloat,bits) do
  e=big"1.704"; gam=BigFloat(Base.MathConstants.eulergamma)
  # Laurent constant of Gamma(z+h)*e^(-z-h) and its z derivatives.
  L(h,p)=begin
   v=BigFloat(factorial(big(h-1)))/e^h; z=sum(inv(BigFloat(j)) for j=1:h-1;init=big(0))-gam-log(e)
   p==0 ? v : p==1 ? v*z : v*(z*z+big(pi)^2/6-sum(inv(BigFloat(j)^2) for j=1:h-1;init=big(0)))
  end
  radial=[L(h,p) for h=1:12,p=0:2]
  allstates=[states(joinpath(get(ENV,"AUDIT_STATE_DIR",joinpath(@__DIR__,"../../claude2/states")),"atoms2_N$(N).txt")) for N in Ns]
  D=zeros(BigFloat,length(Ns)); normcheck=big(0);xq,wq=gauss(nq)
  for (iy,yy) in enumerate(xq),(it,tt) in enumerate(xq)
   y=BigFloat(yy);x=y*BigFloat(tt);hx=sqrt((1+x*x)/2);v=1-2y*y/(1+x*x);Vs=-8/(1-x*x)+1/y
   wt=64/big(pi)*y*(1-x*x)/(1+x*x)^3*v*y*BigFloat(wq[iy])*BigFloat(wq[it]);normcheck+=wt*v
   av=Dict(:one=>big(1),:chi=>chi(x,y),:lR=>log(hx)); xp=[x^j for j=0:12];yp=[y^j for j=0:12]
   for (i,(E,terms)) in enumerate(allstates)
    dp=big(0)
    for (key,g,c) in terms
     n,p,a,b=key;T=c*xp[a+1]*yp[b+1]*av[g]
     if n==0
      dp+=T*(Vs*e*(log(e)+gam-1-log(hx))+E*(gam+log(e)-log(hx)))
     elseif n==1
      dp+=T*(Vs*(-gam-log(e)+log(hx))-E/e)
     else
      dp+=T*(Vs*radial[n-1,p+1]-E*radial[n,p+1])
     end
    end
    D[i]+=wt*dp
   end
  end
  b=-2*(big(pi)-2)/(3big(pi))
  for (i,N) in enumerate(Ns)
   a=-D[i]/4-b/8
   @printf("bits=%d nq=%d N=%d a21=%.15f norm=%.15f\n",bits,nq,N,Float64(a),Float64(normcheck));flush(stdout)
  end
 end
end
if haskey(ENV,"AUDIT_STATE_DIR") && !haskey(ENV,"AUDIT_ALL")
 run(48,256,[6,8,11])
else
 for n in (24,48,96);run(n,192,collect(2:11));end
 run(48,288,[10,11])
end
