# Is the degree-2 log-free part of Ψ, rewritten with ln(L) (L = s, s + r12, R),
# a polynomial in (r1, r2, r12)?  Least-squares test on S^5.
using FockExpansion, LinearAlgebra, Random

const Z = 2.0
const E = -2.903724377034119598311
const A21 = 0.47674787900
const κ = (π - 2) / (3π)

function rows(α, θ, logarg)
    r1, r2 = cos(α / 2), sin(α / 2)  # R = 1
    r12 = sqrt(r1^2 + r2^2 - 2r1 * r2 * cos(θ))
    v = sin(α) * cos(θ)
    L = logarg(r1, r2, r12)
    # ψ20 + ψ21 ln R  with R = 1  ->  write ln R = ln L + ln(R/L): extra ψ21 * ln(1/L)
    target = psi20(α, θ; Z, E, a21 = A21) + (-Z * κ * v) * log(1 / L)
    basis = [r1^2 + r2^2, r1 * r2, r12^2, (r1 + r2) * r12]
    return basis, target
end

function test(name, logarg)
    Random.seed!(1)
    n = 400
    A = zeros(n, 4)
    b = zeros(n)
    for i = 1:n
        α, θ = 0.05 + 3.04rand(), 0.05 + 3.04rand()
        bs, t = rows(α, θ, logarg)
        A[i, :] = bs
        b[i] = t
    end
    x = A \ b
    println(rpad(name, 22), "rms residual = ", norm(A * x - b) / sqrt(n), "   rms target = ", norm(b) / sqrt(n))
end

test("ln s", (r1, r2, r12) -> r1 + r2)
test("ln(s + r12)", (r1, r2, r12) -> r1 + r2 + r12)
test("ln R (Fock)", (r1, r2, r12) -> 1.0)
