using Serialization, LinearAlgebra
include(joinpath(@__DIR__,"../../claude11/extend.jl"))
setprecision(BigFloat,256);T=BigFloat;P=T(pi)
R=deserialize(joinpath(@__DIR__,"x26_chain64.jls"))
O=deserialize(joinpath(@__DIR__,"../../claude11/chain_n64_N48.jls"))
V=inv([T(z)^m for z=0:4,m=0:4])
sector(R,s,kp,m)=sum(V[m+1,z+1]*R.out[(z,s)][kp] for z=0:4)
for s=0:2,kp in ((2,0),(3,0),(4,0),(4,1)),m=1:3
 a=sector(R,s,kp,m);b=sector(O,s,kp,m)
 println("step ",s," ",kp," Z",m," corrected max ",Float64(maximum(abs.(a)))," old-new ",Float64(maximum(abs.(a-b))))
end
ys=R.g.Y
println("corrected double fit error ",maximum(abs.(sector(R,1,(4,0),2).- [1/T(9)-7/T(45)*sin(y+P/2) for y in ys])))
