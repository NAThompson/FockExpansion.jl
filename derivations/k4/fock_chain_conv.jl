# Spectral convergence of the chain up to k = 6: compare grids n = 24, 32, 40 by
# barycentric interpolation at fixed test points. Free harmonics: a₂₁ = 0.4, all others 0.
include("fock_chain.jl")
function bary(xs, f, x)                      # barycentric interpolation, Chebyshev first kind
    n=length(xs); w=[(-1)^j*sin(π*(j+0.5)/n) for j in 0:n-1]
    for j in 1:n; x==xs[j] && return f[j]; end
    s=sum(w[j]/(x-xs[j]) for j in 1:n); sum(w[j]*f[j]/(x-xs[j]) for j in 1:n)/s
end
evalgrid(g, F, X, Y)=bary(g.X, [bary(g.Y, F[i, :], Y) for i in 1:g.n], X)
Z, E = 2.0, -2.9
pts=[(0.3, 2.0), (-1.0, 2.5), (1.2, 1.9), (0.0, π), (-0.5, 1.7)]   # (X, Y) in the rectangle
res=Dict{Int,Any}()
for n in (24, 32, 40)
    g=Grid(n)
    t=@elapsed ψ=fock(g; Z, E, kmax=6, free=Dict((2,1)=>0.4))
    res[n]=Dict(key=>[evalgrid(g, ψ[key], X, Y) for (X, Y) in pts] for key in keys(ψ))
    @printf("n=%d: chain to k=6 in %.2f s\n", n, t)
end
for key in sort(collect(keys(res[40])))
    d1=maximum(abs, res[24][key]-res[40][key]); d2=maximum(abs, res[32][key]-res[40][key])
    @printf("ψ%d%d: |n24-n40| %.1e  |n32-n40| %.1e   value at first point % .10f\n", key..., d1, d2, res[40][key][1])
end
