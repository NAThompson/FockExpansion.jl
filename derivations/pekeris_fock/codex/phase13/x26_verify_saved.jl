using Serialization,LinearAlgebra
include(joinpath(@__DIR__,"../../claude11/extend.jl"))
setprecision(BigFloat,256);T=BigFloat;P=T(pi)
root=joinpath(@__DIR__,"../..")
R=deserialize(joinpath(@__DIR__,"x26_chain80_bounded.jls"));S=deserialize(joinpath(@__DIR__,"x26_chain64.jls"))
O=deserialize(joinpath(root,"claude11/chain_n80_N60.jls"));E=deserialize(joinpath(root,"claude11/ext_side1_n80_N60.jls"))
V=inv([T(z)^m for z=0:4,m=0:4])
sector(R,s,kp,m)=sum(V[m+1,z+1]*R.out[(z,s)][kp] for z=0:4)
function cv(v,x)
 n=length(v);p=T(pi)
 c=[2/T(n)*sum(v[j+1]*cos(p*k*(j+T(1)/2)/n) for j=0:n-1)/(k==0 ? 2 : 1) for k=0:n-1]
 sum(c[k+1]*cos(k*acos(x)) for k=0:n-1)
end
metrics=Pair{String,T}[]
record(k,v)=push!(metrics,k=>T(v))
for kp in ((2,0),(3,0),(4,0),(4,1)),m=0:4
 n=E.g.n
 ext=sum(V[m+1,z+1]*(length(E.results[z].ψ1[kp])>=2 ? E.results[z].ψ1[kp][2][E.e1.iA,:] : zeros(T,n)) for z=0:4)
 a=sector(R,0,kp,m)
 record("corrected_vs_extend_$(kp)_Z$m",maximum(abs.(a-ext)))
 for step=0:2
  a=sector(R,step,kp,m);b=sector(S,step,kp,m)
  record("n64_n80_step$(step)_$(kp)_Z$m",maximum(abs(cv(a,T(x))-cv(b,T(x))) for x in (-.9,-.5,0.,.3,.8)))
 end
end
fit=[1/T(9)-7/T(45)*sin(y+P/2) for y in R.g.Y]
record("corrected_double_fit_max",maximum(abs.(sector(R,1,(4,0),2)-fit)))
record("original_vs_corrected_double",maximum(abs.(sector(R,1,(4,0),2)-sector(O,1,(4,0),2))))
record("psi41_original_minus_twice_corrected",maximum(abs.(sector(O,0,(4,1),2)-2sector(R,0,(4,1),2))))
record("psi40_first_log_original_minus_corrected",maximum(abs.(sector(O,0,(4,0),2)-sector(R,0,(4,0),2))))
record("corrected_double_other_sectors",maximum(maximum(abs.(sector(R,1,kp,m))) for kp in ((2,0),(3,0),(4,0),(4,1)),m=0:4 if (kp,m)!=((4,0),2)))
record("corrected_triple_all",maximum(maximum(abs.(sector(R,2,kp,m))) for kp in ((2,0),(3,0),(4,0),(4,1)),m=0:4))
open(joinpath(@__DIR__,"x26_metrics.tsv"),"w") do io
 for (key,val) in metrics;println(io,key,'\t',val);end
end
for (key,val) in metrics
 startswith(key,"n64_n80") || println(key,": ",val)
end

println("x26 peak_rss_bytes ",Sys.maxrss())
