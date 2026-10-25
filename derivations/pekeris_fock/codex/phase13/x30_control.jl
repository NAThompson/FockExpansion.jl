using DelimitedFiles, LinearAlgebra
setprecision(BigFloat,256); T=BigFloat
D=@__DIR__; R=normpath(joinpath(D,"../.."))
ids=parse.(Int,readlines(joinpath(R,"codex/phase9/x19_fit40_big_ids.txt")))
ls=readlines(joinpath(R,"claude9/spec40.txt"))
pts=[parse.(T,split(ls[i]))[1:2] for i in ids]
Y=readdlm(joinpath(D,"x30_targets80.txt"),T)
Z=readdlm(joinpath(R,"claude11/z1parts.txt"),T)
e3=[-(cos(a/2)+sin(a/2))*(5*(cos(a/2)+sin(a/2))^2-3)/36 for (a,t) in pts]
res=vec(sum(Y,dims=2))+vec(sum(Z,dims=2))-e3
println("fixed coefficient -1 control max_abs=",maximum(abs.(res)))
println("fixed coefficient -1 control held_relative_l2=",norm(res[801:end])/norm(vec(sum(Y[801:end,:],dims=2))))
println("first raw mismatch=",sum(Y[1,:])+sum(Z[1,:])," e3=",e3[1])
# Correct extraction of cubic polynomial's linear coefficient from odd evaluations.
# f(Z)=a1*Z+a3*Z^3 + even terms: (8*(f(1)-f(-1))-(f(2)-f(-2)))/12.
println("closed-form extraction corrected by subtracting e3; no fitted coefficient")
