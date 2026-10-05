# How smooth is ψ₃₀ along slices? Chebyshev coefficient decay on α ∈ (0, π/2) at fixed θ,
# and on θ ∈ (0, π) at fixed α.
using FockExpansion, Printf, LinearAlgebra
P=(Z=2.0, E=-2.9, a21=0.4)
function chebcoef(f, a, b, n)
    x=[cos(π*(j+0.5)/n) for j in 0:n-1]
    y=[f(a+(b-a)*(xi+1)/2) for xi in x]
    [2/n*sum(y[j+1]*cos(k*π*(j+0.5)/n) for j in 0:n-1) for k in 0:n-1]
end
for (name, f, a, b) in (("α-slice θ=1.1", α->psi30(α, 1.1; P..., rtol=1e-14), 1e-3, π/2),
                         ("α-slice θ=0.2", α->psi30(α, 0.2; P..., rtol=1e-14), 1e-3, π/2),
                         ("θ-slice α=0.7", θ->psi30(0.7, θ; P..., rtol=1e-14), 1e-3, π-1e-3),
                         ("θ-slice α=1.5", θ->psi30(1.5, θ; P..., rtol=1e-14), 1e-3, π-1e-3))
    c=chebcoef(f, a, b, 64)
    @printf("%-15s |c_k| at k=8,16,32,48,63: %s\n", name, join([@sprintf("%.1e", abs(c[k+1])) for k in (8,16,32,48,63)], " "))
end
println("through both e-e coalescence points (α = π/2, θ = β):")
for (name, f, a, b) in (("θ ∈ [1e-6, π-1e-6]", θ->psi30(π/2, θ; P..., rtol=1e-14), 1e-6, π-1e-6),)
    c=chebcoef(f, a, b, 96)
    @printf("%-20s |c_k| at k=8,16,32,48,64,95: %s\n", name, join([@sprintf("%.1e", abs(c[k+1])) for k in (8,16,32,48,64,95)], " "))
end
