# Pure solution of (Λ²-32)ψ = h on the angular space, viewed as S³ through
# x = (cos α, sin α cos θ, sin α sin θ cos φ, sin α sin θ sin φ), where Λ² = -4Δ_{S³}.
# The modified Green's function of Δ+8, orthogonal to the degree-2 harmonics, is
#   (Δ+8)G = δ - P₂,   G(γ) = [(γ-π)cos3γ + sin(3γ)/6]/(4π² sin γ),
# so ψ = -(1/4)∫G h dσ is the solution with no Y₄ₗ component, whatever the
# Y₄ₗ content of h. dσ = sin²α sinθ dα dθ dφ.

# r₁₂/R = √(1-sin α cos θ), without cancellation near α = π/2, θ = 0.
xi_stable(α, θ) = sqrt(2sin(π/4-α/2)^2+2sin(α)*sin(θ/2)^2)

function green_k4(omc)              # omc = 1 - cos γ
    γ=2asin(sqrt(min(omc/2, 1.0)))
    s, c=sincos(3γ)
    ((γ-π)*c+s/6)/(4π^2*sin(γ))
end

function graded_points(a, b, toa::Bool, tob::Bool; levels, q = 0.35)
    pts=[a, b]
    toa && append!(pts, a .+ (b-a)/2 .* q .^ (0:levels))
    tob && append!(pts, b .- (b-a)/2 .* q .^ (0:levels))
    sort!(unique!(pts))
end
function split_graded(a, b, x; levels, ends = (true, true))
    a<x<b || return graded_points(a, b, ends...; levels)
    vcat(graded_points(a, x, ends[1], true; levels), graded_points(x, b, true, ends[2]; levels)[2:end])
end
function rule_sum(f, pts, (x, w))
    s=0.0
    for i = 1:(length(pts)-1)
        h=(pts[i+1]-pts[i])/2
        m=(pts[i+1]+pts[i])/2
        for j in eachindex(x)
            s+=w[j]*h*f(m+h*x[j])
        end
    end
    s
end

# Azimuthal integral ∫₀^{2π} G dφ between (α,θ) and (a,t); logarithmic at a = α, t = θ.
# With 1 - cos γ = δ + 2d sin²(φ/2), nearby points (δ < d) give a peak of width √(δ/d)
# at φ = 0. On φ ∈ [0, π/2] the map sin(φ/2) = σ sinh v, σ = √(δ/(2d)), makes
# 1 - cos γ = δ cosh²v and the integrand smooth in v; plain Gauss covers the rest.
const AZIMUTHAL_RULE=gauss(32)
function green_k4_azimuthal(α, θ, a, t, rule = nothing; levels = 0)
    δ=2sin((α-a)/2)^2+2sin(α)*sin(a)*sin((θ-t)/2)^2
    d=sin(α)*sin(θ)*sin(a)*sin(t)
    x, w=AZIMUTHAL_RULE
    s=0.0
    lo=0.0
    if δ<d
        σ=sqrt(δ/(2d))
        vmax=asinh(sin(π/4)/σ)
        for (xi, wi) in zip(x, w)
            v=vmax*(xi+1)/2
            sφ=σ*sinh(v)
            s+=wi*vmax*σ*cosh(v)/sqrt(1-sφ^2)*green_k4(δ*cosh(v)^2)
        end
        lo=π/2
    end
    for (xi, wi) in zip(x, w)
        φ=lo+(π-lo)*(xi+1)/2
        s+=wi*(π-lo)/2*green_k4(δ+2d*sin(φ/2)^2)
    end
    2s
end

# Iterated kernel G₂ = G∘G, the pure solution of (Δ+8)G₂ = G (derivations/k4/g2.py):
#   G₂(γ) = [(γ²-2πγ)/(48π²) + 1/72 + 1/(864π²)] sin 3γ / sin γ, bounded and smooth.
function green2_k4(omc)
    γ=2asin(sqrt(min(omc/2, 1.0)))
    c=cos(γ)
    ((γ^2-2π*γ)/(48π^2)+1/72+1/(864π^2))*(4c^2-1)
end
function green2_k4_azimuthal(α, θ, a, t, rule)
    δ=2sin((α-a)/2)^2+2sin(α)*sin(a)*sin((θ-t)/2)^2
    d=sin(α)*sin(θ)*sin(a)*sin(t)
    2rule_sum(φ->green2_k4(δ+2d*sin(φ/2)^2), (0.0, π/2, Float64(π)), rule)
end

# Pure solution of (Λ²-32)ψ = h, i.e. ψ = -(1/4)G*h. With `h2`, also adds the pure
# solution of (Λ²-32)²χ = h2, i.e. χ = (1/16)G₂*h2.
function solve_k4(h, α, θ; n, levels, h2 = nothing)
    rule=gauss(n)
    A=sort(unique(vcat(split_graded(0.0, π/2, α; levels), split_graded(π/2, Float64(π), α; levels))))
    T=split_graded(0.0, Float64(π), θ; levels)
    function point(a, t)
        v=-h(a, t)*green_k4_azimuthal(α, θ, a, t)/4
        h2===nothing || (v+=h2(a, t)*green2_k4_azimuthal(α, θ, a, t, rule)/16)
        v*sin(t)
    end
    rule_sum(a->sin(a)^2*rule_sum(t->point(a, t), T, rule), A, rule)
end

# Tanh-sinh rule on [a, b]: spectrally accurate for integrands singular at either end.
function de_rule(a, b; h = 0.0625, N = 48)
    xs=Float64[]
    ws=Float64[]
    for j = (-N):N
        τ=j*h
        s=π/2*sinh(τ)
        push!(xs, (a+b)/2+(b-a)/2*tanh(s))
        push!(ws, (b-a)/2*h*π/2*cosh(τ)/cosh(s)^2)
    end
    xs, ws
end
# Same solve for exchange-symmetric sources: the target is folded to α ≤ π/2, sources
# live on α ∈ (0, π/2] with kernel K(x,y) + K(x,ȳ), ȳ the mirror α → π-α, and the domain
# is split at the target. The logarithmic kernel singularity and the Coulomb corners at
# (π/2, 0) and (π/2, π) then sit at rectangle corners, where tanh-sinh converges
# exponentially (h = 0.0625 gives about 1e-14, h = 0.125 about 1e-10).
function solve_k4_symmetric(h, α, θ; h2 = nothing, step = 0.0625)
    α>π/2 && (α=π-α)
    N=ceil(Int, 3/step)
    rule2=gauss(8)
    total=0.0
    for (a0, a1) in ((0.0, α), (α, π/2)), (t0, t1) in ((0.0, θ), (θ, Float64(π)))
        (a1>a0 && t1>t0) || continue
        A=de_rule(a0, a1; h = step, N)
        T=de_rule(t0, t1; h = step, N)
        for (a, wa) in zip(A...), (t, wt) in zip(T...)
            K=green_k4_azimuthal(α, θ, a, t)+green_k4_azimuthal(α, θ, π-a, t)
            v=-h(a, t)*K/4
            if h2!==nothing
                K2=green2_k4_azimuthal(α, θ, a, t, rule2)+green2_k4_azimuthal(α, θ, π-a, t, rule2)
                v+=h2(a, t)*K2/16
            end
            total+=wa*wt*sin(a)^2*sin(t)*v
        end
    end
    total
end

