include("../exterior_factorial/plucker_lib.jl")
using Printf, Random
function main()
setprecision(BigFloat,parse(Int,get(ENV,"AUDIT_PREC","512")))
Random.seed!(13)
b=cuspmodel(2); E=cuspsolve(b).E
ps=plucker_series(b,E,65;prec=512)
println("K d2=",b.K," beta=",Float64(ps.beta))
println("beta minus tr Y1 = ",Float64(ps.beta-tr(riccati_coeffs(b,E,2;prec=512)[2])))
for n in (1,2,20,60)
 r=(ps.C0-ps.sigma*I)*ps.xi[n+1]-((ps.beta-n+1)*I-ps.C1)*ps.xi[n]+(n>=2 ? ps.C2*ps.xi[n-1] : zero(ps.xi[1]))
 println("recurrence relative n=",n," ",Float64(norm(r)/max(norm(ps.xi[n+1]),1)))
end
m=4; A=big.(randn(8,8)); X=big.(randn(8,4)); h=big"1e-60"
println("compound derivative residual ",Float64(norm((wedge(X+h*A*X,ps.S)-wedge(X,ps.S))/h-additive_compound(A,ps.S)*wedge(X,ps.S))))
println("pairing residual ",Float64(abs(plucker_pairing(wedge(X,ps.S),ps.xi[1],ps.S,8)-det([X [Matrix{BigFloat}(I,4,4);slow_solvent(b,E;prec=512)]]))))
# Stream Fock coefficients for two seed vectors, avoiding an N^2 stored tensor.
b=cuspmodel(6); z=cuspsolve(b); m=length(b.K); p=2
println("d6 K=",b.K)
s,f=pencil_roots(b,z.E); println("slow range=",extrema(real(s))," fast range=",extrema(real(f))," max root modulus=",maximum(abs.(vcat(s,f))))
U=big.(b.U); D=big.(b.D); T=big.(b.W)+2big(z.E)*I
N=280; kp=N÷2+3; R=big(10); L=log(R)
seeds=hcat(big.(z.a),big.(randn(m))); seeds[:,2]/=norm(seeds[:,2])
prev=zeros(BigFloat,kp,m,p); prev2=copy(prev); F=zeros(BigFloat,m,p); Q=copy(F); mass=zeros(BigFloat,p)
check=Dict{Int,Matrix{BigFloat}}()
for n=0:N
 c=zeros(BigFloat,kp,m,p); src=copy(c)
 for k=0:n÷2
  n>=1 && (src[k+1,:,:]-=(2(n-1)*U+D)*prev[k+1,:,:]+2(k+1)*U*prev[k+2,:,:])
  n>=2 && (src[k+1,:,:]-=T*prev2[k+1,:,:])
 end
 for i=1:m
  den=n*(n+4)-b.lam[i]
  if den==0
   for k=n÷2:-1:0
    c[k+2,i,:]=(src[k+1,i,:]-(k+1)*(k+2)*c[k+3,i,:])/((2n+4)*(k+1))
   end
   c[1,i,:]=seeds[i,:]
  else
   for k=n÷2:-1:0
    c[k+1,i,:]=(src[k+1,i,:]-(2n+4)*(k+1)*c[k+2,i,:]-(k+1)*(k+2)*c[k+3,i,:])/den
   end
  end
 end
 term=zeros(BigFloat,m,p); dt=copy(term)
 for k=0:n÷2+1
  t=R^n*L^k*c[k+1,:,:]; term+=t
  dt+=(n*t+(k>0 ? k*R^n*L^(k-1)*c[k+1,:,:] : zero(t)))/R
  for j=1:p; mass[j]+=norm(t[:,j]); end
 end
 F+=term; Q+=dt
 if n in (80,100,120,140,160,180,200,220,240,260,280)
  check[n]=vcat(F,Q)
  @printf("N=%d increment F %.3e %.3e D %.3e %.3e\n",n,Float64(norm(term[:,1])),Float64(norm(term[:,2])),Float64(norm(dt[:,1])),Float64(norm(dt[:,2])))
  flush(stdout)
 end
 prev2=prev; prev=c
end
for n in sort(collect(keys(check)))
 @printf("N=%d error to N280 %.3e %.3e\n",n,Float64(norm((check[n]-check[280])[:,1])),Float64(norm((check[n]-check[280])[:,2])))
end
println("sum absolute term norms=",Float64.(mass)," final F norms=",[Float64(norm(F[:,j])) for j=1:p])

end
main()
