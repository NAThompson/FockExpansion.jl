include("goal2_panel_rank.jl")
lines=readlines(joinpath(@__DIR__,"goal2_port_reference.txt")); line=0
for (b,y)=((.3,-.6),(-.27,-.7))
 basepoint[]=(b,y); x=ComplexF64[((k-1)%17-8)/19 for k=1:78]
 for axis=1:2
  global line+=1; vals=parse.(Float64,split(lines[line])); ref=complex.(vals[1:2:end],vals[2:2:end]); actual=connection(b+0im,y+0im,sign(b),axis)*x
  err=maximum(abs.(ref-actual)); @assert err<1e-12; println("connection ",b," ",axis," max error ",err)
 end
 global line+=1; vals=parse.(Float64,split(lines[line])); ref=vals[1:2:end]; actual=real.(outputs(b+0im,y+0im,sign(b))*x)
 err=maximum(abs.(ref-actual)); @assert err<1e-12; println("moments ",b," max error ",err)
end
println("PASS: 312 connection actions and four moment outputs")
