# Phase 2 A: locate the non-polynomial content of R^2 psi20 (and check psi21, psi31, psi42, psi41 Z^1/Z^3/Y4l are polynomial).
using FockExpansion, LinearAlgebra, Random
const Z = 2.0
const E = -2.903724377034119598311
const A21 = 0.47674787900

function point()
    α, θ = 0.05 + 3.04rand(), 0.05 + 3.04rand()
    r1, r2 = cos(α / 2), sin(α / 2)
    r12 = sqrt(1 - sin(α) * cos(θ))
    α, θ, r1, r2, r12
end
monos(r1, r2, r12, d) = [r1^i * r2^j * r12^k for i = 0:d for j = 0:d for k = 0:d if i + j + k == d]

function fitres(f, d; n = 600)
    Random.seed!(2)
    A = zeros(n, length(monos(1, 1, 1, d))); b = zeros(n)
    for i = 1:n
        α, θ, r1, r2, r12 = point()
        A[i, :] = monos(r1, r2, r12, d); b[i] = f(α, θ)
    end
    x = A \ b
    norm(A * x - b) / sqrt(n), norm(b) / sqrt(n)
end


χ(α, θ) = FockExpansion.chi(α, cos(θ))









report(name, f, d) = (r = fitres(f, d); println(rpad(name, 44), "rms residual ", r[1], "   rms ", r[2]))


report("d psi30 / d a21 (difference quotient)", (α, θ) -> psi30(α, θ; Z, E, a21 = 1.0) - psi30(α, θ; Z, E, a21 = 0.0), 3)
report("d psi30 / d E (difference quotient)", (α, θ) -> psi30(α, θ; Z, E = -3.0, a21 = A21) - psi30(α, θ; Z, E = -2.0, a21 = A21), 3)
