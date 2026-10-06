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
# sfock keeps every ψ_kp; at high K keep only what the recurrence needs (k-1, k-2).
function chain(g::SGrid{T}; Z, E, K, out) where {T}
    wV=2sqrt(T(2))*sin.(g.α).*cos.(g.β./2)-4Z*sin.(g.β).*(cos.(g.α./2)+sin.(g.α./2))
    ψ=Dict{Tuple{Int,Int},Matrix{T}}((0, 0)=>ones(T, g.n, g.n))
    z=zeros(T, g.n, g.n)
    get0(k, p)=get(ψ, (k, p), z)
    for k in 1:K
        t=@elapsed begin
            o=SOps(g, k*(k+4))
            H=iseven(k) ? sharmonics(g, k) : Matrix{T}[]
            for p in (k÷2):-1:0
                wh(ψk1)=g.w.*(2(k+2)*(p+1)*ψk1+(p+1)*(p+2)*get0(k, p+2)+2E*get0(k-2, p))-2wV.*get0(k-1, p)
                if iseven(k) && p<k÷2 && haskey(ψ, (k, p+1))
                    G=[sip(g, a, b) for a in H, b in H]
                    s=[sum(g.W./g.w.*h.*wh(get0(k, p+1))) for h in H]
                    cc=-(G\s)./(2(k+2)*(p+1))
                    ψ[(k, p+1)]+=sum(cc[l]*H[l] for l in eachindex(H))
                end
                ψ[(k, p)]=ssolve(g, o, wh(get0(k, p+1)), H)
            end
        end
        # the harmonic part of ψ_k,1 is final only after ψ_k0 is solved, so record now
        rec=Dict(p=>(Float64(maximum(abs, ψ[(k, p)])), [interp(g, ψ[(k, p)], a, θ) for (a, θ) in angs]) for p in 0:k÷2)
        push!(out, (k, t, rec))
        m=maximum(first(v) for v in values(rec))
        @printf("k=%3d  max|ψ_kp| = %.3e  (%.1f s)\n", k, m, t); flush(stdout)
        for kk in collect(keys(ψ)); kk[1]<k-2 && delete!(ψ, kk); end
        k%10==0 && serialize("growth_n$(n)_K$(K).jls", out)
    end
    serialize("growth_n$(n)_K$(K).jls", out)
end
chain(g; Z=2, E=Eex, K, out=Any[])
