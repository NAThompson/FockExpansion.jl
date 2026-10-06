# From the stage-1 run: for each ρ, the largest term max|ψ_kp| ρ^k |log ρ|^p (sets the
# cancellation, hence the digits needed) and the order K at which terms fall below 1e-20.
using DoubleFloats, Serialization, Printf
run=Dict(r[1]=>r for r in deserialize(ARGS[1]))
m(k)=k==0 ? Dict(0=>1.0) : Dict(p=>v[1] for (p, v) in run[k][3])
for ρ in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0)
    terms=[maximum(c*ρ^k*abs(log(ρ))^p for (p, c) in m(k)) for k in 0:120]
    kpk=argmax(terms)-1
    Kc=findlast(>(1e-20*maximum(terms)), terms)-1
    @printf("ρ=%4.1f  largest term %.1e at k=%3d   terms < 1e-20·largest beyond k=%s   last term %.1e\n",
            ρ, maximum(terms), kpk, Kc==120 ? ">120" : string(Kc), terms[end])
end
