include("constrained.jl")
function atomstate(g,x,y)
 g==:chi && return chi4(x,y)
 g==:lR && return log((1+x*x)/2)/2
 atom(g,x,y)
end
function readclaude4(name)
 ls=readlines(joinpath(@__DIR__,"input_states",name));E=parse(BigFloat,split(ls[1],"E_variational = ")[2])
 ts=[(parse.(Int,split(l)[1:4])...,Symbol(split(l)[5]),parse(BigFloat,split(l)[6])) for l in ls[2:end]]
 if get(ENV,"DROP_POLES","0")=="1"
  ts=filter(t->!(t[3]==0 && t[4]==0 && t[5]==:V),ts)
 end
 E,ts
end
function statejet(N)
 E,ts=readclaude4("cmin5_N$(N).txt")
 for (x,y) in ((big".2",big".7"),(big".6",big".8"))
  av=Dict(g=>atomstate(g,x,y) for g in unique(t[5] for t in ts))
  ours=terms((-1,0,0),x,y)
  vals=Dict((n,p)=>sum(c*x^a*y^b*av[g] for (m,r,a,b,g,c) in ts if (m,r)==(n,p);init=big(0)) for (n,p,_,_,_) in ours)
  a=(vals[(2,0)]-ours[3][3])/((1+x*x)/2-y*y)
  err=maximum(abs(vals[(n,p)]-v) for (n,p,v,_,_) in ours if (n,p)!=(2,0))
  println("Claude state N=$N ray=($x,$y) readout=$a other_jet_error=$err")
 end
end
function stateenergy(hinv,N)
 E,ts=readclaude4("cmin5_N$(N).txt");ng=unique(t[5] for t in ts)
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
  av=Dict(g=>atomstate(g,x,y) for g in ng)
  ag=Dict(g=>ForwardDiff.gradient(z->atomstate(g,z[1],z[2]),[x,y]) for g in ng)
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
 @printf("independent state energy constrained min5 N%d hinv=%d E-E0=%.15e E-Eexport=%.4e\n",N,hinv,Float64(H/S-ref),Float64(H/S-E));flush(stdout)
end

if abspath(PROGRAM_FILE)==@__FILE__
for N in parse.(Int,split(get(ENV,"NS","4,8"),","))
 statejet(N)
 for h in parse.(Int,split(get(ENV,"HS","8,12"),","));stateenergy(h,N);end
end

end
