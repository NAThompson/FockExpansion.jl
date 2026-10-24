include("goal2_classical.jl")
lines=readlines(joinpath(@__DIR__,"goal2_classical_reference.txt"));line=0
for (a,p)=((-.4,-1.3),(-.7,-.9))
 global line+=1;x=parse.(Float64,split(lines[line]))
 for axis=1:2
  global line+=1;ref=parse.(Float64,split(lines[line]));actual=classicalC(a,p,0,axis)*x
  err=maximum(abs.(actual-ref));@assert err<2e-13;println("classical derivatives ",a," ",axis," max error ",err)
 end
 global line+=1;ref=parse(Float64,lines[line]);actual=(classicalL(a,p,0)*x)[1]
 err=abs(ref-actual);@assert err<2e-13;println("classical output ",a," error ",err)
end
println("PASS 72 derivative components and two physical classical outputs")
