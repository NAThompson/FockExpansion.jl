# BigFloat evaluator of the Feynman 1-D integral for psi41^(2) (repo src/feynman_k4.jl), tanh-sinh in phi (t = sin^2 phi).
# Returns the raw integral F (= psi41^(2) + c0*Y40 + c2*Y42, c = repo PSI41_Z2_SHIFT, known only in Float64).
using FockExpansion
const FE = FockExpansion
# In-process override (repo file untouched): the repo's Cpow uses (m-1)/m in Float64, ~1e-17 relative error in BigFloat.
FE.eval(:(function Cpow(γ, m)
    s, c = sincos(γ)
    m == 0 && return γ
    m == 1 && return s
    c^(m - 1) * s / m + (m - 1) * Cpow(γ, m - 2) / m
end))

function ts_rule(::Type{T}, h, N) where {T}
    xs = T[]
    ws = T[]
    for j = -N:N
        τ = j * h
        s = T(π) / 2 * sinh(τ)
        push!(xs, (1 + tanh(s)) / 2)          # fraction of [0, π/2]
        push!(ws, h * T(π) / 2 * cosh(τ) / cosh(s)^2 / 2)
    end
    xs, ws
end


# Repo integrand with exact rational constants (the repo's 7/3, 5/3 are Float64 literals: ~1e-17 relative error in BigFloat).
function fi_exact(t, omt, x, σ)
    Flo(v, kind) = FE.Uf(kind, v[1], 1 - v[2], 2v[2] - v[2]^2 - v[1]^2, x)
    Fhi(v, kind) = FE.Uf(kind, σ * (1 - v[1]), v[2], 2v[1] - v[1]^2 - v[2]^2, x)
    F, v0, J = t <= 1 / 2 ? (Flo, [σ * t, t], [1 0; 0 -1]) : (Fhi, [omt, omt], [-σ 0; 0 1])
    J = convert(Matrix{typeof(t)}, J)
    d1 = J * FE.ForwardDiff.gradient(v -> F(v, :log), v0)
    H = J * FE.ForwardDiff.hessian(v -> F(v, :xlog), v0) * J
    -F(v0, :inv) / 3 + 7 * d1[2] / 3 + 5 * H[2, 2] / 3 + σ * H[1, 2]
end

function psi41z2_raw(α::T, θ::T; h = T(1) / 2^5, N = ceil(Int, 4.2 / Float64(h))) where {T}
    x = [cos(α), sin(α) * cos(θ), sin(α) * sin(θ), zero(T)]
    xs, ws = ts_rule(T, h, N)
    s = zero(T)
    for (f, w) in zip(xs, ws)
        f <= 0 || f >= 1 && continue
        φ = T(π) / 2 * f
        t, omt = sin(φ)^2, cos(φ)^2
        (t == 0 || omt == 0) && continue
        s += 2 * w * T(π) / 2 * (fi_exact(t, omt, x, 1) + fi_exact(t, omt, x, -1))
    end
    (T(π) - 2) / (3 * T(π)) / (sqrt(T(2)) * T(π)) * s
end

if abspath(PROGRAM_FILE) == @__FILE__
    for (α, θ) in ((0.7, 1.1), (1.3, 2.4), (0.2, 0.3))
        r = FE.psi41_z2(α, θ) + FE.PSI41_Z2_SHIFT[1] * FE.harmonic40(α, θ) + FE.PSI41_Z2_SHIFT[2] * FE.harmonic42(α, θ)
        for bits in (128, 256)
            setprecision(BigFloat, bits) do
                a, t = BigFloat(α), BigFloat(θ)
                v1 = psi41z2_raw(a, t; h = BigFloat(1) / 16)
                v2 = psi41z2_raw(a, t; h = BigFloat(1) / 32)
                println((α, θ), " bits=$bits  F(h/32)=", Float64(v2), "  |h16-h32|=", Float64(abs(v1 - v2)), "  vs repo ", Float64(v2) - r)
            end
        end
    end
end
