include("constrained.jl")
using FockExpansion
mx=big(0)
for (x,y) in ((big".2",big".7"),(big".6",big".8"),(big".05",big".3"))
 q=(1+x*x)/2;al=2atan((1-x)/(1+x));theta=acos((q-y*y)/((1-x*x)/2))
 p1,p2,L=jets(x,y)
 ref=(sqrt(q)*FockExpansion.psi10(al,theta;Z=big(2)),q*FockExpansion.psi20(al,theta;Z=big(2),E=E0,a21=big(0)),q*FockExpansion.psi21(al,theta;Z=big(2)))
 println("ours=",(p1,p2,L)," ref=",ref," chi=",chi4(x,y)," refchi=",FockExpansion.chi(al,cos(theta)))
 er=maximum(abs.([p1,p2,L].-collect(ref)));global mx=max(mx,er)
 println("source angular jet comparison x=$x y=$y max_error=$er")
end
@assert mx<big"1e-45"
raycheck()
println("PASS seed convention and radial expansion checks")
