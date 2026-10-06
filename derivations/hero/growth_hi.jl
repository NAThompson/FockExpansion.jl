# Stage 1 of the hero run: the bare chain (free constants zero) to high order in Double64.
# Records max|ψ_kp| and values at fixed off-grid angles, to measure decay and grid resolution.
# usage: julia --project=. growth_hi.jl n K
using DoubleFloats, Serialization
include("../k4/fock_spectral.jl")
const T=Double64
n=parse(Int, ARGS[1]); K=parse(Int, ARGS[2])
const Eex=T(big"-2.903724377034119598311159245194404446696905")
angs=[(T(π)*(i-T(1)/2)/8, T(π)*(j-T(1)/2)/8) for i in 1:8 for j in 1:8]
g=SGrid(T, n)
out=Any[]
function record(k, ψ, t)
    rec=Dict(p=>(Float64(maximum(abs, ψ[(k, p)])), [interp(g, ψ[(k, p)], a, θ) for (a, θ) in angs]) for p in 0:k÷2)
    push!(out, (k, t, rec))
    @printf("k=%3d  max|ψ_kp| = %.3e  (%.1f s)\n", k, maximum(first(v) for v in values(rec)), t); flush(stdout)
    k%10==0 && serialize("growth_n$(n)_K$(K).jls", out)
end
sfock(g; Z=2, E=Eex, kmax=K, keep=2, onk=record)
serialize("growth_n$(n)_K$(K).jls", out)
