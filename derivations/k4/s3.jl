# Quadrature over the (α, θ) square with the Fock measure sin²α sinθ dα dθ.
# Panels are graded geometrically toward α = 0, π/2, π and θ = 0, where the
# Coulomb factors 1/sin(α/2), 1/cos(α/2) and 1/ξ are singular.
using QuadGK
function graded(a, b, toa::Bool, tob::Bool; levels=12, q=0.35)
    pts=[a, b]
    if toa; append!(pts, a .+ (b-a)/2 .* q.^(0:levels)); end
    if tob; append!(pts, b .- (b-a)/2 .* q.^(0:levels)); end
    sort!(unique!(pts))
end
const GL=Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()
rule(n)=get!(()->gauss(n), GL, n)
function panelsum(f, pts, n)
    x, w=rule(n); s=0.0
    for i in 1:length(pts)-1
        a, b=pts[i], pts[i+1]; h=(b-a)/2; m=(a+b)/2
        for j in eachindex(x); s+=w[j]*h*f(m+h*x[j]); end
    end
    s
end
apts(levels)=sort(unique(vcat(graded(0, π/2, true, true; levels), graded(π/2, π, true, true; levels))))
tpts(levels)=graded(0, π, true, false; levels)
const APTS=apts(12)
const TPTS=tpts(12)
sphint(f; n=24, levels=12)=panelsum(α->sin(α)^2*panelsum(θ->f(α, θ)*sin(θ), tpts(levels), n), apts(levels), n)
