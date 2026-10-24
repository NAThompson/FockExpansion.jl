include("../exterior_factorial/exterior_lib.jl")
using Printf
function d2check()
 b=cuspmodel(2); E=cuspsolve(b).E
 N=800; Ys=riccati_coeffs(b,E,N+2;prec=3072)
 for w in (0.224,0.3)
  B=factorial_coeffs(Ys,w,N)
  for R in (0.75,1.0,2.0)
   ref=cuspexterior(b,E,R;rmax=30.,order=30,fraction=.3)
   ps=factorial_partial_sums(Ys[1],B,w,R)
   @printf("d2 w=%.3f R=%.2f errors ",w,R)
   for k in (100,200,300,400,600,800); @printf("K%d=%.6e ",k,norm(ps[k+1]-ref)); end
   println()
  end
 end
end
function d6check()
 b=cuspmodel(6); z=cuspsolve(b); E=z.E
 Y1=cuspexterior(b,E,9.;rmax=30.,order=30,fraction=.3)
 ref=cuspexterior(b,E,1.;rmax=30.,order=30,fraction=.3)
 ys=taylor_connect(b,E,Y1,9.,1.,600)
 c=cuspcoeff(b,E,44); F,D,_=evalseries(c,1.); M=D/F
 ev=eigen(M-ref); v=ev.vectors[:,argmin(real(ev.values))]
 for n in (320,420,460,520,600)
  @printf("d6 exact-panel start Taylor n=%d maperr=%.8e null-action=%.8e mismatch=%.8e\n",n,norm(ys[n]-ref),norm((ys[n]-ref)*v),minimum(real(eigvals(M-ys[n]))))
 end
 @printf("d6 panels mismatch=%.8e\n",minimum(real(ev.values)))
end
function endpointcheck()
 for d in (2,4,6)
  b=cuspmodel(d); E=cuspsolve(b).E
  Ys=riccati_coeffs(b,E,60;prec=512)
  for R in (12.,24.)
   a=asymptotic_partial_sums(Ys,R); yr=cuspexterior(b,E,R;rmax=30.,order=30,fraction=.3)
   @printf("d%d R%.0f 2term-vs-high=%.6e asym40-vs-panel30=%.6e default-vs-panel30=%.6e\n",d,R,norm(a[2]-a[41]),norm(a[41]-yr),norm(cuspexterior(b,E,R)-yr))
  end
 end
end
d2check(); flush(stdout); d6check(); flush(stdout); endpointcheck()
