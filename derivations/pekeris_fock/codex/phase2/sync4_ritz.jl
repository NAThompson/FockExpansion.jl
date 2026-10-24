# Read-only reuse of Claude's assembly, with independent driver and diagnostics.
setprecision(BigFloat,parse(Int,get(ENV,"AUDIT_BITS","384")))
include(joinpath(@__DIR__,"../../claude2/atoms_ritz.jl"))
const ae=big"1.704"
hinv=parse(Int,get(ENV,"AUDIT_HINV","24"));nd=nodes(inv(BigFloat(hinv)),5hinv)
at=Dict(:one=>(ones(BigFloat,length(nd.x)),zeros(BigFloat,length(nd.x)),zeros(BigFloat,length(nd.x))))
at[:chi]=atom_table(chi_xy,nd);at[:lR]=atom_table(lR_xy,nd)
println("bits=",precision(BigFloat)," hinv=",hinv," nodes=",length(nd.x));flush(stdout)
M=moments(nd,at,28,28);R=Radial(ae)
println("moment error=",M.tab[(:one,:v,:one,:v)][3,4]/2-big(1)/21);flush(stdout)
for N in (6,8,11)
 B=build(N,:atoms2;ε=ae);H,S=matrices(B,M,R,ae);En,c=ground(H,S)
 C=cholesky(Symmetric(S));res=norm(H*c-En*S*c)/(norm(H)*norm(c)+abs(En)*norm(S)*norm(c))
 @printf("N=%d size=%d delta=%.16e scaled_res=%.4e cond=%.4e\n",N,length(B),Float64(En-E0),Float64(res),Float64(condest(S)));flush(stdout)
 outdir=joinpath(@__DIR__,"audit$(hinv)");mkpath(outdir)
 open(joinpath(outdir,"atoms2_N$(N).txt"),"w") do io
  println(io,"# eps = 1.704 E_variational = ",En)
  acc=Dict{Tuple{Int,Int,Int,Int,Symbol},BigFloat}()
  for (bf,cb) in zip(B,c), t in bf, ((a,b,g),cy) in t.Y
   key=(t.n,t.p,a,b,g);acc[key]=get(acc,key,big(0))+cb*cy/c[1]
  end
  for (key,cf) in sort(collect(acc);by=first);println(io,join(key," ")," ",cf);end
 end
end
