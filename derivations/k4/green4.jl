# Pure solution of (Λ² - 32)ψ = h on S³ (k = 4, resonant), via the modified
# Green's function of Δ + 8 orthogonal to the degree-2 harmonics:
#   (Δ+8)G = δ - P₂,   G(γ) = [(γ-π)cos3γ + sin3γ/6]/(4π² sinγ),
# so ψ(x) = -(1/4)∫ G(x·y) h(y) dσ(y) and ψ ⊥ Y₄ₗ automatically.
include("s3.jl")
function Gpure(omc)               # omc = 1 - cos γ, accurate near γ = 0
    γ=2asin(sqrt(min(omc/2, 1.0)))
    s, c=sincos(3γ)
    ((γ-π)*c+s/6)/(4π^2*sin(γ))
end
# Azimuthal integral ∫₀^{2π} G dφ between (α,θ) and (α',θ').
function kernel(α, θ, a, t; levels=14, n=8)
    δ=2sin((α-a)/2)^2+2sin(α)*sin(a)*sin((θ-t)/2)^2   # 1 - cos γ at φ = 0
    d=sin(α)*sin(θ)*sin(a)*sin(t)
    pts=graded(0.0, Float64(π), true, false; levels, q=0.3)
    2panelsum(φ->Gpure(δ+2d*sin(φ/2)^2), pts, n)
end
function split_graded(a, b, x; levels=12, q=0.35, ends=(true, true))
    if a<x<b
        vcat(graded(a, x, ends[1], true; levels, q), graded(x, b, true, ends[2]; levels, q)[2:end])
    else
        graded(a, b, ends...; levels, q)
    end
end
function solve4(h, α, θ; n=12, levels=12)
    A=sort(unique(vcat(split_graded(0.0, π/2, α; levels), split_graded(π/2, Float64(π), α; levels))))
    T=split_graded(0.0, Float64(π), θ; levels, ends=(true, true))
    -panelsum(a->sin(a)^2*panelsum(t->h(a, t)*sin(t)*kernel(α, θ, a, t), T, n), A, n)/4
end
