# Green's-function solve with double-exponential (tanh-sinh) tensor rules.
# Target folded to α ≤ π/2 (sources are exchange symmetric), source domain α ∈ (0, π/2]
# with kernel K(x,y) + K(x,ȳ), ȳ the mirror α → π-α; the domain is split at the target,
# so the logarithmic kernel singularity and the Coulomb corners (π/2, 0), (π/2, π) all
# sit at rectangle corners, where tanh-sinh converges spectrally.
using FockExpansion, Printf
const F=FockExpansion
function de_rule(a, b; h=0.125, N=24)
    xs=Float64[]; ws=Float64[]
    for j in -N:N
        τ=j*h; s=π/2*sinh(τ)
        x=tanh(s); w=h*π/2*cosh(τ)/cosh(s)^2
        push!(xs, (a+b)/2+(b-a)/2*x); push!(ws, (b-a)/2*w)
    end
    xs, ws
end
function solve_de(h, α, θ; h2=nothing, kw...)
    α>π/2 && (α=π-α)
    total=0.0
    for (a0, a1) in ((0.0, α), (α, π/2)), (t0, t1) in ((0.0, θ), (θ, Float64(π)))
        a1>a0 && t1>t0 || continue
        A=de_rule(a0, a1; kw...); T=de_rule(t0, t1; kw...)
        for (a, wa) in zip(A...), (t, wt) in zip(T...)
            K=F.green_k4_azimuthal(α, θ, a, t)+F.green_k4_azimuthal(α, θ, π-a, t)
            v=-h(a, t)*K/4
            if h2!==nothing
                K2=F.green2_k4_azimuthal(α, θ, a, t, F.gauss(8))+F.green2_k4_azimuthal(α, θ, π-a, t, F.gauss(8))
                v+=h2(a, t)*K2/16
            end
            total+=wa*wt*sin(a)^2*sin(t)*v
        end
    end
    total
end
