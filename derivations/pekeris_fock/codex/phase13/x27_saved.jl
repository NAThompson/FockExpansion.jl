using Serialization,LinearAlgebra
setprecision(BigFloat,384)
root=joinpath(@__DIR__,"../..")
R=deserialize(joinpath(root,"claude11/lad_B1.jls"))
S=deserialize(joinpath(root,"claude11/lad_B1_low.jls"))
function ev(c,x)
 b1=zero(x);b2=zero(x)
 for k=length(c)-1:-1:1
  b1,b2=2*x*b1-b2+c[k+1],b1
 end
 x*b1-b2+c[1]
end
println("keys ",sort(collect(keys(R)),by=string))
c=R[(0,(4,0),"A")]
println("length ",length(c)," precision_first ",precision(c[1]))
println("tail_last5 ",[string(x) for x in c[end-4:end]])
for ty in [-3//4,-1//4,0//1,1//4,3//4]
 x=BigFloat(ty);u=cos(5*BigFloat(pi)/4+BigFloat(pi)*x/4)
 println("ty=",ty," u=",u," B_code=",ev(c,x)," A3_code=",ev(S[(0,(3,0),"A")],x))
end
include(joinpath(root,"claude9/fock_spectral_copy.jl"))
T=deserialize(joinpath(root,"claude11/lad_fixed_0-1101010.jls"))
println("saved chain steps ",T.steps," grid ",T.g.n)
for st=0:T.steps
 println("step ",st," max A4 ",maximum(abs.(T.out[(0,st)][(4,0)]))," max A3 ",maximum(abs.(T.out[(0,st)][(3,0)])))
end
println("maxrss_bytes ",Sys.maxrss())
