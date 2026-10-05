include("fock_chain.jl")
using FockExpansion
Z, E, a21 = 2.0, -2.9, 0.4
g=Grid(24)
# a₂₁ is the coefficient of sin α cos θ = x₂, the l = 1 harmonic at k = 2 (C₀⁽²⁾ = 1, sinα P₁ = x₂)
ψ=fock(g; Z, E, kmax=4, free=Dict((2, 1)=>a21))
nodes=[αθ(g.X[i], g.Y[j]) for i in 1:g.n, j in 1:g.n]
pkg=Dict((1,0)=>(a,t)->psi10(a,t; Z), (2,1)=>(a,t)->psi21(a,t; Z), (2,0)=>(a,t)->psi20(a,t; Z, E, a21),
         (3,1)=>(a,t)->psi31(a,t; Z), (3,0)=>(a,t)->psi30(a,t; Z, E, a21, rtol=1e-14),
         (4,2)=>(a,t)->psi42(a,t; Z), (4,1)=>(a,t)->psi41(a,t; Z, E, a21))
for key in sort(collect(keys(pkg)))
    err=maximum(abs(ψ[key][i, j]-pkg[key](nodes[i, j]...)) for i in 1:g.n, j in 1:g.n)
    @printf("ψ%d%d: max |chain - package| over the grid = %.1e\n", key..., err)
end
idx=((5,9), (12,12), (18,4))
err40=maximum(abs(ψ[(4,0)][i,j]-psi40(nodes[i,j]...; Z, E, a21)) for (i,j) in idx)
@printf("ψ40: max |chain - package| at 3 nodes = %.1e\n", err40)
