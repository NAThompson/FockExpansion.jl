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
function green_k4_azimuthal(α, θ, a, t, rule; levels)
    δ=2sin((α-a)/2)^2+2sin(α)*sin(a)*sin((θ-t)/2)^2
    d=sin(α)*sin(θ)*sin(a)*sin(t)
    2rule_sum(φ->green_k4(δ+2d*sin(φ/2)^2), graded_points(0.0, Float64(π), true, false; levels = levels+6, q = 0.3), rule)
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
        v=-h(a, t)*green_k4_azimuthal(α, θ, a, t, rule; levels)/4
        h2===nothing || (v+=h2(a, t)*green2_k4_azimuthal(α, θ, a, t, rule)/16)
        v*sin(t)
    end
    rule_sum(a->sin(a)^2*rule_sum(t->point(a, t), T, rule), A, rule)
end
