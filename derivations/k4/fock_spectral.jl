# Precision-generic spectral solution of the Fock hierarchy (see fock_chain.jl for the
# method). Eigenpairs of the separated 1D operators are computed in Float64 and refined by
# inverse iteration in the working precision, so only LU (generic in LinearAlgebra) is needed.
using LinearAlgebra, Printf
struct SGrid{T}
    n::Int; X::Vector{T}; Y::Vector{T}; α::Matrix{T}; β::Matrix{T}; w::Matrix{T}; W::Matrix{T}; D::Matrix{T}
end
function SGrid(::Type{T}, n) where {T}
    P=T(π)
    x=[cos(P*(j+T(1)/2)/n) for j in 0:n-1]
    bw=[(-1)^j*sin(P*(j+T(1)/2)/n) for j in 0:n-1]
    D=zeros(T, n, n)
    for i in 1:n, j in 1:n
        i==j || (D[i, j]=bw[j]/bw[i]/(x[i]-x[j]))
    end
    for i in 1:n; D[i, i]=-sum(D[i, :]); end
    X=P/2 .* x; Y=P .+ P/2 .* x
    α=[(X[i]+Y[j])/2 for i in 1:n, j in 1:n]; β=[(Y[j]-X[i])/2 for i in 1:n, j in 1:n]
    w=2sin.(α).*sin.(β)
    θs=[P*(j+T(1)/2)/n for j in 0:n-1]
    q=[2/T(n)*(1-2sum(cos(2k*t)/(4k^2-1) for k in 1:n÷2)) for t in θs].*(P/2)
    W=[(P/2)*w[i, j]*q[i]*q[j] for i in 1:n, j in 1:n]
    SGrid{T}(n, X, Y, α, β, w, W, D)
end
# Newton refinement of the whole eigendecomposition M = VΛV⁻¹ (eigenvalues distinct):
# with A = V⁻¹MV, the correction V ← V(I+E), E_ij = A_ij/(λ_j-λ_i), converges quadratically.
# Each step is a few n³ operations in the working precision.
function refine_eigen(M::Matrix{T}) where {T}
    e=eigen(Float64.(M))
    V=T.(real(e.vectors))
    λ=T.(real(e.values))
    tol=eps(T)*opnorm(M, 1)
    prev=T(Inf)
    for it in 1:20
        A=V\(M*V)
        λ=diag(A)
        off=maximum(abs(A[i, j]) for i in axes(A, 1), j in axes(A, 2) if i!=j)
        # stop at roundoff: below tolerance, or no longer converging quadratically
        (off<tol || off>prev/4) && break
        prev=off
        E=[i==j ? zero(T) : A[i, j]/(λ[j]-λ[i]) for i in axes(A, 1), j in axes(A, 2)]
        V+=V*E
        V./=sqrt.(sum(abs2, V; dims=1))
        it==20 && @warn "refine_eigen: off-diagonal $(Float64(off)) after 20 steps"
    end
    λ, V
end
struct SOps{T}; VX::Matrix{T}; λX::Vector{T}; iVX::Matrix{T}; VY::Matrix{T}; λY::Vector{T}; iVY::Matrix{T}; end
function SOps(g::SGrid{T}, c) where {T}
    DX=(2/T(π))*g.D
    MX=Diagonal(cos.(g.X))*DX^2-Diagonal(sin.(g.X))*DX+(T(c)/16)*Diagonal(cos.(g.X))
    MY=Diagonal(cos.(g.Y))*DX^2-Diagonal(sin.(g.Y))*DX+(T(c)/16)*Diagonal(cos.(g.Y))
    λX, VX=refine_eigen(Matrix(MX)); λY, VY=refine_eigen(Matrix(MY))
    SOps{T}(VX, λX, inv(VX), VY, λY, inv(VY))
end
function sharmonics(g::SGrid{T}, k) where {T}
    x1=cos.(g.α); x2=cos.(g.β); s2=1 .- x1.^2
    H=Matrix{T}[]
    for l in 0:k÷2
        m=k÷2-l
        C0=ones(T, size(x1)); C1=2(l+1)*x1; Cm=m==0 ? C0 : C1
        for j in 2:m
            Cn=(2(j+l)*x1.*C1-(j+2l)*C0)/j; C0, C1=C1, Cn; Cm=Cn
        end
        # sin^l α P_l(cos θ) = Q_l(x2, s2): (j+1)Q_{j+1} = (2j+1)x2 Q_j - j s2 Q_{j-1}
        Q0=ones(T, size(x1)); Pl=l==0 ? Q0 : x2
        Q1=x2
        for j in 1:l-1
            Qn=((2j+1)*x2.*Q1-j*s2.*Q0)/(j+1); Q0, Q1=Q1, Qn; Pl=Qn
        end
        push!(H, Cm.*Pl)
    end
    H
end
sip(g, f, h)=sum(g.W.*f.*h)
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
# The harmonics as columns of an n²×m matrix with a factored Gram matrix, built once per k:
# projections are then one matrix-vector product instead of m² grid sums.
struct HBasis{T}; H::Vector{Matrix{T}}; Hm::Matrix{T}; WHm::Matrix{T}; G::LU{T,Matrix{T},Vector{Int}}; end
function HBasis(g::SGrid{T}, k) where {T}
    H=iseven(k) ? sharmonics(g, k) : Matrix{T}[]
    Hm=isempty(H) ? zeros(T, g.n^2, 0) : reduce(hcat, vec.(H))
    WHm=vec(g.W).*Hm
    HBasis{T}(H, Hm, WHm, lu(Hm'*WHm))
end
Base.isempty(hb::HBasis)=isempty(hb.H)
# coefficients c with ⟨H_l, f⟩ = Σ_m G_lm c_m, for the inner product with weight W·u
hcoef(hb::HBasis, f, u=nothing)=hb.G\(hb.WHm'*(u===nothing ? vec(f) : vec(f).*vec(u)))
hcombo(hb::HBasis, c)=reshape(hb.Hm*c, size(hb.H[1]))
function ssolve(g::SGrid{T}, o::SOps{T}, wh, hb::HBasis{T}) where {T}
    Rt=o.iVX*(-wh./16)*o.iVY'
    res=resonant_mask(o.λX, o.λY, length(hb.H))
    F=[res[i, j] ? zero(T) : Rt[i, j]/(o.λX[i]-o.λY[j]) for i in 1:g.n, j in 1:g.n]
    ψ=o.VX*F*o.VY'
    isempty(hb) && return ψ
    ψ-hcombo(hb, hcoef(hb, ψ))
end
# `onk(k, ψ, seconds)` is called after each order; with `keep = 2` only the two latest orders
# (all the recurrence needs) are kept, for long runs.
function sfock(g::SGrid{T}; Z, E, kmax, free=Dict{Tuple{Int,Int},Any}(), keep=typemax(Int), onk=nothing) where {T}
    Z, E=T(Z), T(E)
    wV=2sqrt(T(2))*sin.(g.α).*cos.(g.β./2)-4Z*sin.(g.β).*(cos.(g.α./2)+sin.(g.α./2))
    ψ=Dict{Tuple{Int,Int},Matrix{T}}((0, 0)=>ones(T, g.n, g.n))
    get0(k, p)=get(ψ, (k, p), zeros(T, g.n, g.n))
    for k in 1:kmax
        t0=time()
        o=SOps(g, k*(k+4))
        hb=HBasis(g, k)
        for p in (k÷2):-1:0
            wh(ψk1)=g.w.*(2(k+2)*(p+1)*ψk1+(p+1)*(p+2)*get0(k, p+2)+2E*get0(k-2, p))-2wV.*get0(k-1, p)
            if iseven(k) && p<k÷2 && haskey(ψ, (k, p+1))
                cc=-hcoef(hb, wh(get0(k, p+1)), 1 ./ g.w)./(2(k+2)*(p+1))
                ψ[(k, p+1)]+=hcombo(hb, cc)
            end
            ψkp=ssolve(g, o, wh(get0(k, p+1)), hb)
            iseven(k) && p==0 && (ψkp+=sum(T(get(free, (k, l), 0))*hb.H[l+1] for l in 0:k÷2))
            ψ[(k, p)]=ψkp
        end
        onk===nothing || onk(k, ψ, time()-t0)
        for kk in collect(keys(ψ)); kk[1]<=k-keep && delete!(ψ, kk); end
    end
    ψ
end
# Barycentric interpolation of a grid function at (α, θ).
function interp(g::SGrid{T}, F, α, θ) where {T}
    P=T(π)
    β=acos(sin(α)*cos(θ))
    tx, ty=2(α-β)/P, 2(α+β-P)/P
    x=[cos(P*(j+T(1)/2)/g.n) for j in 0:g.n-1]
    bw=[(-1)^j*sin(P*(j+T(1)/2)/g.n) for j in 0:g.n-1]
    ax=bw./(tx.-x); ay=bw./(ty.-x)
    (ax'*F*ay)/(sum(ax)*sum(ay))
end
