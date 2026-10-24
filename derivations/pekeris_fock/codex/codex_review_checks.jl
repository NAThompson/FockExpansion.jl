include(joinpath(@__DIR__,"../julia/galerkin_log.jl"))
for bits in (256,512)
 setprecision(BigFloat,bits) do
  r=solve(6,:fixedlog;bits)
  e=big"1.704"; B=r.B; M=Moments(60)
  P0=[apply_P(f,e,0,2,1) for f in B]
  Q=[mergewith(-,apply_P(f,e,1,2,1),P0[j]) for (j,f) in enumerate(B)]
  H=[-inner(M,f,g) for f in B,g in P0]; S=[inner(M,f,g) for f in B,g in Q]
  cholesky(Hermitian(S))
  residual=norm(H*r.x-r.E*S*r.x)/(norm(H)*norm(r.x)+abs(r.E)*norm(S)*norm(r.x))
  @printf("N=6 bits=%d E=%.30f n=%d asym=%.3e residual=%.3e Gram_cholesky=PASS\n",bits,r.E,r.n,r.asym,residual); flush(stdout)
 end
end
r=solve(13,:fixedlog;bits=768)
@printf("N=13 bits=768 E=%.35f n=%d asym=%.3e\n",r.E,r.n,r.asym)
