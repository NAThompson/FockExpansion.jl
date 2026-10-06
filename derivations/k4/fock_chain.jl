# The whole Fock hierarchy from ψ₀₀ = 1 by separated spectral solves (no closed forms).
#   [Λ² - k(k+4)]ψ_kp = h_kp,  h_kp = 2(k+2)(p+1)ψ_k,p+1 + (p+1)(p+2)ψ_k,p+2 - 2Vψ_k-1,p + 2Eψ_k-2,p.
# Coordinates X = α-β, Y = α+β separate Λ² (sylvester30.jl). Even k is resonant on the
# harmonics Y_kl (l = 0..k/2); their coefficients in ψ_kp (p ≥ 1) are fixed by solvability
# of the ψ_k,p-1 equation, and in ψ_k0 they are free inputs (a₂₁ at k = 2, a₄₀, a₄₂ at k = 4).
using LinearAlgebra, Printf
include("sylvester30.jl")
function fejer(n)
    θ=[π*(j+0.5)/n for j in 0:n-1]
    [2/n*(1-2sum(cos(2k*θj)/(4k^2-1) for k in 1:n÷2)) for θj in θ]
end
struct Grid
    n::Int; X::Vector{Float64}; Y::Vector{Float64}
    α::Matrix{Float64}; β::Matrix{Float64}; w::Matrix{Float64}; W::Matrix{Float64}
end
function Grid(n)
    x, _=chebdiff(n); X=π/2 .* x; Y=π .+ π/2 .* x
    α=[(X[i]+Y[j])/2 for i in 1:n, j in 1:n]; β=[(Y[j]-X[i])/2 for i in 1:n, j in 1:n]
    w=2sin.(α).*sin.(β)                                   # cos X - cos Y
    q=fejer(n).*(π/2)
    W=[(π/2)*w[i, j]*q[i]*q[j] for i in 1:n, j in 1:n]   # S³ measure / (2π)… up to a constant
    Grid(n, X, Y, α, β, w, W)
end
# S-state harmonics of degree k (k even) as functions of x₁ = cos α, x₂ = cos β:
# zonal Gegenbauer in x₁ times sin^l α P_l(cos θ), sin α cos θ = x₂.
function harmonics(g::Grid, k)
    x1=cos.(g.α); x2=cos.(g.β); s2=1 .- x1.^2
    H=Matrix{Float64}[]
    for l in 0:k÷2
        m=k÷2-l
        # C_m^{(l+1)}(x₁) by recurrence
        C0=ones(size(x1)); C1=2(l+1)*x1
        Cm=m==0 ? C0 : C1
        for j in 2:m
            Cm, C1, C0=(2(j+l)*x1.*C1-(j+2l)*C0)/j, (2(j+l)*x1.*C1-(j+2l)*C0)/j, C1
        end
        # sin^l α P_l(cos θ) as a polynomial in x₂ and sin²α:  Σ_k c_k x₂^{l-2k} (sin²α)^k
        Pl=zeros(size(x1))
        for kk in 0:l÷2
            coef=(-1)^kk*binomial(l, kk)*binomial(2l-2kk, l)/2^l
            Pl.+=coef.*x2.^(l-2kk).*s2.^kk
        end
        push!(H, Cm.*Pl)
    end
    H
end
struct Ops
    VX; λX; iVX; VY; λY; iVY
end
function Ops(g::Grid, c)
    _, _, MX, MY=separated_ops(g.n, c)
    eX=eigen(Matrix(MX)); eY=eigen(Matrix(MY))
    VX=real(eX.vectors); VY=real(eY.vectors)
    Ops(VX, real(eX.values), inv(VX), VY, real(eY.values), inv(VY))
end
ip(g::Grid, f, h)=sum(g.W.*f.*h)
# Resonant pairs λX_i = λY_j of the Sylvester operator: exactly the nres = k÷2+1 S-state
# harmonics for even k and none for odd k. Take the nres smallest gaps and check that they
# are separated from the rest, so a near-crossing cannot be mistaken for a harmonic.
function resonant_mask(λX, λY, nres)
    d=[abs(x-y) for x in λX, y in λY]
    p=sortperm(vec(d))
    nres==0 || d[p[nres]]<1e-6*d[p[nres+1]] || error("resonant gaps not separated: $(d[p[nres]]) vs $(d[p[nres+1]])")
    d[p[nres+1]]>1e-2 || error("near-resonant non-harmonic pair, gap $(d[p[nres+1]])")
    m=falses(size(d)); m[p[1:nres]].=true
    m
end
# pure solution of (Λ² - c)ψ = h given wh = (cos X - cos Y)h on the grid
function solve(g::Grid, o::Ops, wh, H)
    Rt=o.iVX*(-wh./16)*o.iVY'
    res=resonant_mask(o.λX, o.λY, length(H))
    F=[res[i, j] ? 0.0 : Rt[i, j]/(o.λX[i]-o.λY[j]) for i in 1:g.n, j in 1:g.n]
    ψ=o.VX*F*o.VY'
    isempty(H) && return ψ
    G=[ip(g, a, b) for a in H, b in H]
    ψ-sum(c*h for (c, h) in zip(G\[ip(g, h, ψ) for h in H], H))
end
"""ψ[(k,p)] on the grid for k ≤ kmax. `free[(k,l)]` = coefficient of the l-th harmonic in ψ_k0."""
function fock(g::Grid; Z, E, kmax, free=Dict{Tuple{Int,Int},Float64}())
    wV=2sqrt(2)*sin.(g.α).*cos.(g.β./2)-4Z*sin.(g.β).*(cos.(g.α./2)+sin.(g.α./2))   # (cosX-cosY)·V
    ψ=Dict{Tuple{Int,Int},Matrix{Float64}}((0, 0)=>ones(g.n, g.n))
    get0(k, p)=get(ψ, (k, p), zeros(g.n, g.n))
    for k in 1:kmax
        c=k*(k+4); o=Ops(g, c)
        H=iseven(k) ? harmonics(g, k) : Matrix{Float64}[]
        for p in (k÷2):-1:0
            # (cos X - cos Y)·h_kp, without the harmonic part of ψ_k,p+1 still to be fixed
            wh(ψk1)=g.w.*(2(k+2)*(p+1)*ψk1+(p+1)*(p+2)*get0(k, p+2)+2E*get0(k-2, p))-2wV.*get0(k-1, p)
            ψkp=solve(g, o, wh(get0(k, p+1)), H)
            if iseven(k) && p<k÷2 && haskey(ψ, (k, p+1))
                # fix the harmonic part of ψ_k,p+1 by solvability of this (p) equation:
                # ⟨Y, h_kp⟩ = 0, h_kp linear in ψ_k,p+1 with coefficient 2(k+2)(p+1)
                G=[ip(g, a, b) for a in H, b in H]
                s=[sum(g.W./g.w.*h.*wh(get0(k, p+1))) for h in H]
                cc=-(G\s)./(2(k+2)*(p+1))
                ψ[(k, p+1)]+=sum(cc[l]*H[l] for l in eachindex(H))
                ψkp=solve(g, o, wh(ψ[(k, p+1)]), H)
            end
            if iseven(k) && p==0
                ψkp+=sum(get(free, (k, l), 0.0)*H[l+1] for l in 0:k÷2)
            end
            ψ[(k, p)]=ψkp
        end
    end
    ψ
end
