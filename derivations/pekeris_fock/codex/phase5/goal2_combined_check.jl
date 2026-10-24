include("goal2_combined.jl")
lines=readlines(joinpath(@__DIR__,"goal2_combined_reference.txt"));line=0
function parsevec(line)
 v=parse.(Float64,split(line));complex.(v[1:2:end],v[2:2:end])
end
setprecision(160) do
 for center=((-.4,-1.3),(-.6,-1.1))
  a,p=big.(center);global line+=1;x=parsevec(lines[line]);global line+=1;z=parsevec(lines[line])
  err=maximum(abs.(x-physical_projection(a,p)*z));@assert err<2e-13;println("actual endpoint quotient ",center," error ",Float64(err))
  for axis=1:2
   global line+=1;ref=parsevec(lines[line]);C,l=combined(complex(a),complex(p),center,axis)
   actual=C*x;err=maximum(abs.(actual[331:332]-ref));@assert err<2e-12;println("old atom pullback ",center," ",axis," error ",Float64(err))
  end
 end
end
println("PASS actual endpoint quotients and old-atom pullbacks")
