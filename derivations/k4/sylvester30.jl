# Spectral solution of the Fock recurrence in separating coordinates.
# X = α - β, Y = α + β (β = angle to e₂, cos β = sin α cos θ) map the physical region to
# the rectangle X ∈ [-π/2, π/2], Y ∈ [π/2, 3π/2], with the coalescences at the corners, and
#   Δ_{S³} = 4/(cos X - cos Y) [∂_X(cos X ∂_X) - ∂_Y(cos Y ∂_Y)].
# So (Λ² - c)ψ = h, Λ² = -4Δ, becomes the separated equation
#   [∂_X cosX ∂_X + (c/16)cos X]ψ - [∂_Y cosY ∂_Y + (c/16)cos Y]ψ = -(cos X - cos Y) h/16,
# a Sylvester equation on a tensor Chebyshev grid. Validated here on ψ₃₀ (c = 21).
using FockExpansion, LinearAlgebra, Printf
const F=FockExpansion
# Chebyshev points of the first kind on [-1,1] and their differentiation matrix
function chebdiff(n)
    x=[cos(π*(j+0.5)/n) for j in 0:n-1]
    w=[(-1)^j*sin(π*(j+0.5)/n) for j in 0:n-1]          # barycentric weights
    D=zeros(n, n)
    for i in 1:n, j in 1:n
        i==j || (D[i, j]=w[j]/w[i]/(x[i]-x[j]))
    end
    for i in 1:n; D[i, i]=-sum(D[i, :]); end
    x, D
end
function separated_ops(n, c)
    x, D=chebdiff(n)
    X=π/2 .* x; Y=π .+ π/2 .* x
    DX=(2/π)*D
    MX=Diagonal(cos.(X))*DX^2-Diagonal(sin.(X))*DX+(c/16)*Diagonal(cos.(X))
    MY=Diagonal(cos.(Y))*DX^2-Diagonal(sin.(Y))*DX+(c/16)*Diagonal(cos.(Y))
    X, Y, MX, MY
end
αθ(X, Y)=(α=(X+Y)/2; β=(Y-X)/2; (α, acos(clamp(cos(β)/sin(α), -1.0, 1.0))))
Z, E, a21 = 2.0, -2.9, 0.4
# r = -(cos X - cos Y) h₃₀/16, with h₃₀ = 10ψ₃₁ - 2Vψ₂₀ + 2Eψ₁₀; cos X - cos Y = 2 sin α sin β
# cancels the Coulomb singularities: 2 sinα sinβ/ξ = 2√2 sinα cos(β/2), sinα V₁ = 2(cos(α/2)+sin(α/2)).
function rhs30(X, Y)
    α=(X+Y)/2; β=(Y-X)/2; _, θ=αθ(X, Y)
    w=2sin(α)*sin(β)
    wV=2sqrt(2)*sin(α)*cos(β/2)-Z*sin(β)*2*2*(cos(α/2)+sin(α/2))     # w·V
    h_noV=10psi31(α, θ; Z)+2E*psi10(α, θ; Z)
    -(w*h_noV-2wV*psi20(α, θ; Z, E, a21))/16
end
for n in (12, 16, 24, 32)
    X, Y, MX, MY=separated_ops(n, 21.0)
    R=[rhs30(X[i], Y[j]) for i in 1:n, j in 1:n]
    Ψ=sylvester(Matrix(MX), -Matrix(MY)', -R)
    err=maximum(abs(Ψ[i, j]-psi30(αθ(X[i], Y[j])...; Z, E, a21, rtol=1e-14)) for i in 1:n, j in 1:n)
    @printf("n=%2d: max |spectral - closed form| on the grid = %.1e\n", n, err)
end
