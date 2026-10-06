# Independent check of the chain: apply Λ² in the original coordinates (α, θ) to the
# barycentric interpolant (exact derivatives by ForwardDiff) and compare with h_kp
# at off-grid points, for every coefficient up to k = 6.
include("fock_chain.jl")
using ForwardDiff
function bary(xs, f, x)
    n=length(xs); w=[(-1)^j*sin(π*(j+0.5)/n) for j in 0:n-1]
    s=sum(w[j]/(x-xs[j]) for j in 1:n); sum(w[j]*f[j]/(x-xs[j]) for j in 1:n)/s
end
evalgrid(g, F, X, Y)=bary(g.X, [bary(g.Y, F[i, :], Y) for i in 1:g.n], X)
# value at (α, θ): β = acos(sin α cos θ), X = α-β, Y = α+β
val(g, F, α, θ)=(β=acos(sin(α)*cos(θ)); evalgrid(g, F, α-β, α+β))
function Λ2(g, F, α, θ)
    f(v)=val(g, F, v[1], v[2])
    H=ForwardDiff.hessian(f, [α, θ]); ∇=ForwardDiff.gradient(f, [α, θ])
    -4*(H[1,1]+2cot(α)*∇[1]+(H[2,2]+cot(θ)*∇[2])/sin(α)^2)
end
Z, E = 2.0, -2.9
g=Grid(32)
ψ=fock(g; Z, E, kmax=6, free=Dict((2,1)=>0.4, (4,0)=>0.1, (4,2)=>-0.2))
V(α, θ)=1/sqrt(1-sin(α)*cos(θ))-Z*(1/sin(α/2)+1/cos(α/2))
for (α, θ) in ((0.7, 1.1), (2.2, 2.4), (1.4, 0.35))
    worst=0.0
    for (k, p) in sort(collect(keys(ψ)))
        k==0 && continue
        get0(kk, pp)=haskey(ψ, (kk, pp)) ? val(g, ψ[(kk, pp)], α, θ) : 0.0
        h=2(k+2)*(p+1)*get0(k, p+1)+(p+1)*(p+2)*get0(k, p+2)-2V(α, θ)*get0(k-1, p)+2E*get0(k-2, p)
        r=Λ2(g, ψ[(k, p)], α, θ)-k*(k+4)*val(g, ψ[(k, p)], α, θ)-h
        worst=max(worst, abs(r))
    end
    @printf("(α,θ)=(%.1f,%.2f): max over all ψ_kp (k ≤ 6) of |(Λ²-k(k+4))ψ_kp - h_kp| = %.1e\n", α, θ, worst)
end
