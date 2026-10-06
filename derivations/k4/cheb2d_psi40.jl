# 2D Chebyshev interpolation of ψ₄₀ (tanh-sinh evaluator at ~1e-10) in (α, s), with β = arccos(sin α cos θ) = π/2 + α s,
# α ∈ (0, π/2], s ∈ [-1, 1]. Tests spectral convergence over the whole physical region.
using FockExpansion, Printf, Random
P=(Z=2.0, E=-2.9, a21=0.4)
θof(α, s)=acos(clamp(cos(π/2+α*s)/sin(α), -1.0, 1.0))
f(α, s)=psi40(α, θof(α, s); P..., step=0.125)
cheb_nodes(n)=[cos(π*(j+0.5)/n) for j in 0:n-1]
function fit(n)
    xa=cheb_nodes(n); xs=cheb_nodes(n)
    A(x)=π/4*(x+1)                         # α ∈ (0, π/2)
    V=[f(A(xa[i]), xs[j]) for i in 1:n, j in 1:n]
    T(k, x)=cos(k*acos(x))
    C=[ (2/n)^2*sum(V[i,j]*T(k,xa[i])*T(l,xs[j]) for i in 1:n, j in 1:n) for k in 0:n-1, l in 0:n-1]
    C[1,:]./=2; C[:,1]./=2
    C
end
evalc(C, x, y)=(n=size(C,1); sum(C[k+1,l+1]*cos(k*acos(x))*cos(l*acos(y)) for k in 0:n-1, l in 0:n-1))
Random.seed!(1)
test=[(rand(), 2rand()-1) for _ in 1:40]
for n in (8, 12, 16)
    C=fit(n)
    err=maximum(abs(evalc(C, 2a-1, s)-f(π/4*(2a), s)) for (a, s) in test)
    @printf("n=%2d: max interpolation error over 40 random points %.1e;  |C| tail %.1e\n", n, err, maximum(abs, C[end, :]) + maximum(abs, C[:, end]))
end
