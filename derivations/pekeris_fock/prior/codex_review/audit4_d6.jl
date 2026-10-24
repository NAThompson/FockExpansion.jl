include("../exterior_factorial/exterior_lib.jl")
include_string(Main, join(split(first(split(read("../exterior_factorial/twoseries.jl", String), "base = Dict(")), "\n")[4:end], "\n"))
b=cuspmodel(6); E=big(cuspsolve(b).E); h=big"1e-20"
p=(Rs=10.0,K=240,omega=0.284,nfock=210,prec=384,rprec=1024)
f0=match(b,E;p...).f; s0=(match(b,E+h;p...).f-f0)/h
q=merge(p,(Rs=11.0,nfock=235,)); f1=match(b,E;q...).f; s1=(match(b,E+h;q...).f-f1)/h
println("base slope=",s0," radius slope=",s1," ratio=",s1/s0)
println("estimated R residual correction from supplied dE=",big"-4.550097227542718e-36"*(s0/s1-1))
