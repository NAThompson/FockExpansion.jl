include("x19_atoms.jl")
data=[parse.(Float64,split(l)) for l in readlines(joinpath(X19DIR,"../../claude9/spec40.txt"))]
ids=findall(r->0.25<r[1]<π-0.25 && 0.25<r[2]<π-0.25,data)
rules=[let (z,w)=gaussleg(n);((z.+1)./2,w./2) end for n in (64,128,256,512)]
worst=(0.0,0,0)
for i in ids
    g=geometry(data[i][1:2]...);a=panels(g...,rules[2]...);b=panels(g...,rules[3]...)
    er,j=findmax(abs.((a-b)./(1 .+abs.(b))))
    er>worst[1] && (global worst=(er,i,j))
end
println("worst discrepancy, input row, panel-state index: ",worst)
g=geometry(data[worst[2]][1:2]...);println("geometry ",g)
for (n,r) in zip((64,128,256,512),rules)
    v=panels(g...,r...);println("n=",n," value=",v[worst[3]])
end
