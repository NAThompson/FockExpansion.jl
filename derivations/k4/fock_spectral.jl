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
function refine_eigen(M::Matrix{T}) where {T}
    e=eigen(Float64.(M))
    λs=T.(real(e.values)); V=T.(real(e.vectors))
    # Fixed shift at the Float64 eigenvalue: one LU per pair, and each step gains about
    # -log10(1e-14/gap) digits. A Rayleigh shift would make M-λI exactly singular.
    iters=ceil(Int, precision(T)*log10(2)/8)+1
    for k in eachindex(λs)
        v=V[:, k]
        F=lu(M-(λs[k]*(1+T(2)^-40)+T(2)^-40)*I)
        for _ in 1:iters
            y=F\v; v=y/norm(y)
        end
        V[:, k]=v; λs[k]=dot(v, M*v)
    end
    λs, V
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
        Pl=zeros(T, size(x1))
        for kk in 0:l÷2
            Pl.+=T((-1)^kk*binomial(l, kk)*binomial(2l-2kk, l))/2^l .* x2.^(l-2kk).*s2.^kk
        end
        push!(H, Cm.*Pl)
    end
    H
end
sip(g, f, h)=sum(g.W.*f.*h)
function ssolve(g::SGrid{T}, o::SOps{T}, wh, H) where {T}
    Rt=o.iVX*(-wh./16)*o.iVY'
    F=[abs(o.λX[i]-o.λY[j])<1e-6 ? zero(T) : Rt[i, j]/(o.λX[i]-o.λY[j]) for i in 1:g.n, j in 1:g.n]
    ψ=o.VX*F*o.VY'
    isempty(H) && return ψ
    G=[sip(g, a, b) for a in H, b in H]
    ψ-sum(c*h for (c, h) in zip(G\[sip(g, h, ψ) for h in H], H))
end
function sfock(g::SGrid{T}; Z, E, kmax, free=Dict{Tuple{Int,Int},Any}()) where {T}
    Z, E=T(Z), T(E)
    wV=2sqrt(T(2))*sin.(g.α).*cos.(g.β./2)-4Z*sin.(g.β).*(cos.(g.α./2)+sin.(g.α./2))
    ψ=Dict{Tuple{Int,Int},Matrix{T}}((0, 0)=>ones(T, g.n, g.n))
    get0(k, p)=get(ψ, (k, p), zeros(T, g.n, g.n))
    for k in 1:kmax
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
            ψkp=ssolve(g, o, wh(get0(k, p+1)), H)
            iseven(k) && p==0 && (ψkp+=sum(T(get(free, (k, l), 0))*H[l+1] for l in 0:k÷2))
            ψ[(k, p)]=ψkp
        end
    end
    ψ
end
