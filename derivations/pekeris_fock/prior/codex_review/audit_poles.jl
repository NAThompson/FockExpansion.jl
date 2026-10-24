include("../exterior_factorial/exterior_lib.jl")
using Printf
function polecheck(d,p0,h)
 b=cuspmodel(d); E=cuspsolve(b).E; Ys=riccati_coeffs(b,E,45;prec=256)
 m=length(b.K); U=b.U; D=b.D; T=b.W+2E*I; L=Diagonal(b.lam)
 r=14+im*imag(p0); Y=sum(ComplexF64.(Ys[n+1])/r^n for n=0:40)
 f(r,Y)=-(Y*Y+(5/r*I+2U)*Y-L/r^2+D/r+T)
 while real(r)>4+1e-10
  step=-min(h,real(r)-4)
  k1=f(r,Y); k2=f(r+step/2,Y+step/2*k1); k3=f(r+step/2,Y+step/2*k2); k4=f(r+step,Y+step*k3)
  Y+=step/6*(k1+2k2+2k3+k4); r+=step
 end
 eye=Matrix{ComplexF64}(I,m,m); N=120; B=[zeros(ComplexF64,m,m) for _=0:N]; B[1]=eye; B[2]=r*Y
 for n=0:N-2
  rhs=((n+1)*(2n+5)*eye+2r*(n+1)*U)*B[n+2]
  rhs+=(n*(n+4)*eye+4r*n*U-L+r*D+r^2*T)*B[n+1]
  n>=1 && (rhs+=(2r*(n-1)*U+r*D+2r^2*T)*B[n])
  n>=2 && (rhs+=r^2*T*B[n-1])
  B[n+3]=-rhs/((n+1)*(n+2))
 end
 p=p0; sv=0.
 for _=1:12
  x=(p-r)/r; P=copy(B[end]); Q=zero(P)
  for n=N-1:-1:0; Q=Q*x+P; P=P*x+B[n+1]; end
  F=svd(P); v=F.V[:,end]; w=F.U[:,end]; sv=F.S[end]/F.S[1]
  dp=-dot(w,P*v)/dot(w,(Q/r)*v); p+=dp
  abs(dp)<1e-13 && break
 end
 @printf("d%d h=%.4f pole=%.12f%+.12fi sigma-ratio=%.3e\n",d,h,real(p),imag(p),sv)
end
for (d,p) in ((2,.636+1.85im),(4,1.404+4.16im),(6,.56+3.47im))
 for h in (.004,.002); polecheck(d,p,h); flush(stdout); end
end
