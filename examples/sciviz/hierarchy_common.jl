# Shared code for the hierarchy_*.jl figures: the spectral chain solver of
# derivations/k4/fock_spectral.jl (separating coordinates X = α - β, Y = α + β), plus
# barycentric interpolation of its grid functions to (α, θ) and S-state harmonics.
include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "..", "..", "derivations", "k4", "fock_spectral.jl"))
using ForwardDiff, Printf
const latexstring = CairoMakie.Makie.LaTeXStrings.latexstring

# Categorical colours (fixed order) for line plots.
const SERIES =
    ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6d4bc4", "#7a7a7a"]
const DIVERGING = Reverse(:RdBu)

"""All ψ_kp for k ≤ kmax at the helium ground-state inputs. a₂₁ is the only free
constant set; the free harmonic parts of ψ_k0 for even k ≥ 4 are zero."""
function he_chain(::Type{T} = Float64; n = 64, kmax = 24) where {T}
    g = SGrid(T, n)
    ψ = sfock(g; Z = T(HE.Z), E = T(HE.E), kmax, free = Dict((2, 1) => T(HE.a21)))
    g, ψ
end

keys_upto(kmax) = [(k, p) for k = 0:kmax for p = 0:(k÷2)]

"""Chebyshev (first kind) nodes and barycentric weights used by `SGrid`."""
function cheb_nodes(::Type{T}, n) where {T}
    P = T(π)
    x = [cos(P * (j + T(1) / 2) / n) for j = 0:(n-1)]
    bw = [(-1)^j * sin(P * (j + T(1) / 2) / n) for j = 0:(n-1)]
    x, bw
end

"""Barycentric weights of the point (α, θ) in the X and Y directions."""
function bary_weights(x, bw, α, θ)
    P = oftype(one(α) * π, π)
    β = acos(sin(α) * cos(θ))
    tx = 2(α - β) / P
    ty = 2(α + β - P) / P
    ax = bw ./ (tx .- x)
    ay = bw ./ (ty .- x)
    ax ./ sum(ax), ay ./ sum(ay)
end

"""Value of the grid function F at (α, θ). Generic, so ForwardDiff can differentiate it."""
function interp(x, bw, F, α, θ)
    ax, ay = bary_weights(x, bw, α, θ)
    ax' * F * ay
end

"""Interpolate every coefficient in ψ to the (αs, θs) grid. Returns Dict of [iα, iθ]."""
function interp_grid(g, ψ, ks, αs, θs)
    x, bw = cheb_nodes(Float64, g.n)
    Fs = [Float64.(ψ[key]) for key in ks]
    out = [Matrix{Float64}(undef, length(αs), length(θs)) for _ in ks]
    Threads.@threads for j in eachindex(θs)
        for i in eachindex(αs)
            ax, ay = bary_weights(x, bw, αs[i], θs[j])
            for (m, F) in enumerate(Fs)
                out[m][i, j] = ax' * F * ay
            end
        end
    end
    Dict(zip(ks, out))
end

"""The angular potential V(α, θ) (README)."""
Vpot(α, θ; Z = HE.Z) = 1 / sqrt(1 - sin(α) * cos(θ)) - Z * (1 / cos(α / 2) + 1 / sin(α / 2))

"""Λ²f at (α, θ) by ForwardDiff, for f(α, θ)."""
function Λ2(f, α, θ)
    h(v) = f(v[1], v[2])
    v = [α, θ]
    ∇ = ForwardDiff.gradient(h, v)
    H = ForwardDiff.hessian(h, v)
    -4 * (H[1, 1] + 2cot(α) * ∇[1] + (H[2, 2] + cot(θ) * ∇[2]) / sin(α)^2)
end

"""Right-hand side h_kp of the recurrence, from a value lookup `val(k, p)` (0 if absent)."""
recurrence_rhs(val, k, p, V, E) =
    2(k + 2) * (p + 1) * val(k, p + 1) + (p + 1) * (p + 2) * val(k, p + 2) -
    2V * val(k - 1, p) + 2E * val(k - 2, p)

"""Singlet S-state hyperspherical harmonics of Λ²-degree K = 2n (n = S³ degree):
C_m^(l+1)(cos α) sin^l α P_l(cos θ), m + l = n, m even, on the chain grid (cos β =
sin α cos θ). Built by three-term recurrences (no binomials), so large K is fine.
Returns (l, Y) pairs."""
function singlet_harmonics(g, K)
    n = K ÷ 2
    x1 = cos.(g.α)
    x2 = cos.(g.β)
    s2 = 1 .- x1 .^ 2
    out = Tuple{Int,Matrix{Float64}}[]
    Qm, Q = zero(x1), one.(x1)          # sin^l α P_l(cos θ): Q_{l-1}, Q_l
    for l = 0:n
        if l > 0
            Qm, Q = Q, ((2l - 1) .* x2 .* Q .- (l - 1) .* s2 .* Qm) ./ l
        end
        m = n - l
        isodd(m) && continue
        λ = l + 1
        C0, C1 = one.(x1), 2λ .* x1
        Cm = m == 0 ? C0 : C1
        for j = 2:m
            Cm = (2(j + λ - 1) .* x1 .* C1 .- (j + 2λ - 2) .* C0) ./ j
            C0, C1 = C1, Cm
        end
        push!(out, (l, Cm .* Q))
    end
    out
end

"""Mean square of f over S³ with the chain quadrature."""
rms(g, f) = sqrt(sum(g.W .* f .^ 2) / sum(g.W))

label_kp(k, p) = latexstring("\\psi_{$k,$p}")
