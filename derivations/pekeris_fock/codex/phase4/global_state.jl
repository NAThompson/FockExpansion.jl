include("review_state.jl")
function globalstate(N,nq)
 E,ts=readclaude4("cmin5_N$(N).txt");z,w=gauss(nq);I=big(0);nc=big(0)
 rad=Dict((m,p)=>moment(m,p,ee) for m=0:12 for p=0:2)
 for i=1:nq,j=1:nq
  y=z[i];x=y*z[j];q=(1+x*x)/2;v=1-y*y/q;V=-8/(1-x*x)+1/y
  wt=64/big(pi)*y*y*(1-x*x)/(1+x*x)^3*v*w[i]*w[j];nc+=wt*v
  av=Dict(g=>atomstate(g,x,y) for g in unique(t[5] for t in ts));dp=big(0)
  for (n,p,a,b,g,c) in ts
   dp+=c*x^a*y^b*av[g]*(n==0 ? V*ee*(log(ee)+gamma-1-log(q)/2)+E*(gamma+log(ee)-log(q)/2) : n==1 ? V*(-gamma-log(ee)+log(q)/2)-E/ee : V*rad[(n-2,p)]-E*rad[(n-1,p)])
  end
  I+=wt*dp
 end
 println("GLOBAL Claude N=$N nq=$nq a21=$(-I/4+2kap/8) projection_norm=$nc");flush(stdout)
end
for N in parse.(Int,split(get(ENV,"NS","4,8,10,11"),","))
 statejet(N)
 for q in (24,48);globalstate(N,q);end
end
