# Read-only audit of Claude arrays. No solver or analysis script imported.
using Serialization, LinearAlgebra, Printf
include(joinpath(@__DIR__, "../../claude11/extend.jl"))
setprecision(BigFloat,256)
T=BigFloat; P=T(pi)
R=deserialize(joinpath(@__DIR__,"../../claude11/ext_side1_n80_N60.jls"))
C=deserialize(joinpath(@__DIR__,"../../claude11/chain_n80_N60.jls"))
D=deserialize(joinpath(@__DIR__,"../../claude11/chain_n64_N48.jls"))
V=inv([T(z)^m for z=0:4,m=0:4])
sector(C,s,kp,m)=sum(V[m+1,z+1]*C.out[(z,s)][kp] for z=0:4)
n=R.g.n
ys=R.g.Y
qy=[2/T(n)*(1-2sum(cos(2k*t)/(4k^2-1) for k=1:n÷2))*P/2 for t in [P*(j+T(1)/2)/n for j=0:n-1]]
a1=sum(V[2,z+1]*(length(R.results[z].ψ1[(2,0)])>=2 ? R.results[z].ψ1[(2,0)][2][R.e1.iA,:] : zeros(T,n)) for z=0:4)
println("psi20 Z1 code A1 projection onto cos(Y/2): ",sum(qy.*a1.*cos.(ys./2))/(P/2))
println("psi20 Z1 code A1 projection onto sin(Y/2): ",sum(qy.*a1.*sin.(ys./2))/(P/2))
for m=0:4
 b=sector(C,1,(4,0),m)
 target=m==2 ? [1/T(9)-7/T(45)*sin(y+P/2) for y in ys] : zeros(T,n)
 println("double variation sector ",m," error ",maximum(abs.(b-target)))
end
for kp in ((2,0),(3,0),(4,0),(4,1)),m=0:4
 a=sector(C,2,kp,m)
 println("triple ",kp," Z",m,": ",maximum(abs.(a)))
end
println("Saved-data check completed")
cb=sector(C,0,(2,0),1)
println("chain vs extend psi20 Z1 far-log difference: ",maximum(abs.(cb-a1)))
println("chain psi20 Z1 sin projection: ",sum(qy.*cb.*sin.(ys./2))/(P/2))
# Check already-final harmonic boundary copied in chain step 0.
g=SGrid(T,16);ph=sfock(g;Z=one(T),E=zero(T),kmax=2)
println("physical psi21 max = ",maximum(abs.(ph[(2,1)])))
