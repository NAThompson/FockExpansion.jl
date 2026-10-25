using LinearAlgebra, Printf
BLAS.set_num_threads(1)
include(joinpath(@__DIR__, "../../claude9/fock_spectral_copy.jl"))
function weights(n,x)
 b=[(-1.0)^j*sin(pi*(j+.5)/n) for j=0:n-1]
 function at(t)
  k=findfirst(z->abs(z-t)<1e-14,x)
  if k!==nothing
   a=zeros(n); a[k]=1; return a
  end
  a=b./(t.-x); a/sum(a)
 end
 at
end
function simp(f,a,b,n=32768)
 d=(b-a)/n
 d/3*(f(a)+f(b)+sum((isodd(k) ? 4 : 2)*f(a+k*d) for k=1:n-1))
end
function f1(a,x,y)
 # F1(a;1/2,1/2;a+1/2;x,y), a=r+1/2
 total=0.0; cx=1.0
 for m=0:100
  m>0 && (cx*=((m-.5)/m)*x)
  cy=1.0
  for n=0:100
   n>0 && (cy*=((n-.5)/n)*y)
   ratio=1.0
   for j=0:m+n-1; ratio*=(a+j)/(a+.5+j); end
   total+=cx*cy*ratio
  end
 end
 total
end
A=.2; D=.15; c=A+D
p=[4c^3-3c, -24c*c*D+6D,48c*D^2,-32D^3]
x=-2D/(1-c); y=2D/(1+c)
rat=1.; total=0.
for r=0:3
 global rat,total
 r>0 && (rat*=((r-.5)/r))
 total+=p[r+1]*rat*f1(r+.5,x,y)
end
appell=2pi/sqrt(1-c*c)*total
angular=simp(t->(4(A+D*cos(t))^3-3(A+D*cos(t)))/sqrt(1-(A+D*cos(t))^2),0.,2pi)
@printf("elliptic angular %.17g Appell %.17g error %.3g\n",angular,appell,abs(angular-appell))
@assert abs(angular-appell)<1e-12
z0=.37; gamma=complex(acos(z0)); initial=gamma
for k=1:2048
 z=1+(z0-1)*exp(2pi*im*k/2048)
 g=acos(z)
 candidates=[s*g+2pi*j for s in (-1,1), j=-2:2]
 global gamma=candidates[argmin(abs.(candidates.-gamma))]
end
G(g)=((g-pi)*cos(3g)+sin(3g)/6)/(4pi^2*sin(g))
delta=G(gamma)-G(initial); expected=(4z0^3-3z0)/(2pi*sqrt(1-z0^2))
@printf("DeltaG loop %s expected %.17g error %.3g\n",string(delta),expected,abs(delta-expected))
@assert abs(delta-expected)<1e-13
for n in (24,40,56)
 g=SGrid(Float64,n); o=SOps(g,32)
 xx=cos.(pi.*((0:n-1).+.5)./n); at=weights(n,xx)
 b=at(1.); ex=at((3.4-pi)*2/pi); ey=at((3.0-pi)*2/pi); ee=at((2.8-pi)*2/pi)
 V=o.VY; vb=vec(b'*V); vx=vec(ex'*V); vy=vec(ey'*V)
 # Discrete product operator on nodal source vectors (quadrature is in inverse V).
 P=V*Diagonal(vx./vb)*o.iVY
 err=0.
 for h in (sin.(g.Y),1 .+sqrt(2)*cos.(g.Y),1 .-sqrt(2)*cos.(g.Y))
  target=(dot(ex,h)/dot(b,h)).*h
  err=max(err,maximum(abs.(P*h-target)))
 end
 # Point kernel obtained by integrating source with cardinal quadrature weights.
 q=[simp(t->begin a=at((t-pi)*2/pi); a[k] end,pi/2,3pi/2,2048) for k=1:n]
 kv=dot(ey,P*(ee./q))
 endpoint=norm(V*Diagonal(vb./vb)*o.iVY-I,Inf)
 @printf("SOps n=%d weak_exact_modes_error=%.3g endpoint_identity=%.3g point_discrete=%.17g\n",n,err,endpoint,kv)
 @assert err<1e-7 && endpoint<1e-9
end
println("CHECKS PASS; point_discrete is diagnostic, NOT a certified continuum value")
println("peak_rss_bytes=",Sys.maxrss())
