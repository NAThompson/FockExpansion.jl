using LinearAlgebra, Printf
include(joinpath(@__DIR__,"../julia/pekeris33.jl"))
for cutoff in (12,16,20)
 A,B=pencil(cutoff,Float64); ev=eigen(A,-B)
 good=[i for i in eachindex(ev.values) if isfinite(ev.values[i]) && abs(imag(ev.values[i]))<1e-8]
 j=good[argmax(real.(ev.values[good]))]; e=real(ev.values[j]); v=real.(ev.vectors[:,j]); idx=index_set(cutoff)
 v ./= sqrt(sum((t[1]==t[2] ? 1 : 2)*v[i]^2 for (i,t) in enumerate(idx)))
 @printf("cutoff=%d eps=%.15f relative_residual=%.3e\n",cutoff,e,norm((A+e*B)*v)/(norm(A+e*B)*norm(v)))
 for N=0:cutoff
  z=sqrt(sum((t[1]==t[2] ? 1 : 2)*v[i]^2 for (i,t) in enumerate(idx) if sum(t)==N))
  @printf("%d %.10e\n",N,z)
 end
end
println("Exact 1D projection coefficients of exp((1/2-beta)*x); beta<0 has a divergent Laguerre expansion")
for beta in (0.3944,-0.1)
 t=beta+0.5; ratio=(t-1)/t
 @printf("beta=%g ratio=%g |a_10/a_0|=%g\n",beta,ratio,abs(ratio)^10)
end
