# Stage 4(b): Rayleigh–Ritz in the Pekeris space at fixed eps, with and without a log(sigma) sector.
# Functions are exp(-sigma/2) f, f = sum c * u^a v^b w^c sigma^-q log(sigma)^p.
# Matrix entries <f_i, P f_j> use the polynomial operator P = -W'(H - E) generated from sympy
# (P(E) = P0 + E Q, Q = W' > 0), and exact moments
#   ∫ u^a v^b w^c σ^-q log^p σ e^-σ = a! b! c! Γ(N+3-q)/Γ(N+3) * {1, ψ, ψ^2 + ψ'}(N+3-q).
using LinearAlgebra, Printf
include("generated.jl")

const Key = NTuple{5,Int}  # (a, b, c, p, q)
Term = Dict{Key,BigFloat}

addto!(d::Term, k::Key, x) = (d[k] = get(d, k, zero(BigFloat)) + x; d)

function deriv(f::Term, i::Int)
    g = Term()
    for (k, c) in f
        a = collect(k)
        if a[i] > 0
            b = copy(a); b[i] -= 1
            addto!(g, Tuple(b), c * a[i])
        end
        k[5] > 0 && addto!(g, (k[1], k[2], k[3], k[4], k[5] + 1), -c * k[5])
        k[4] > 0 && addto!(g, (k[1], k[2], k[3], k[4] - 1, k[5] + 1), c * k[4])
    end
    g
end

function apply_P(f::Term, epsv, E, Z, lam)
    out = Term()
    cache = Dict{NTuple{3,Int},Term}((0, 0, 0) => f)
    getd(d) = get!(cache, d) do
        # differentiate step by step
        base = d[1] > 0 ? (d[1] - 1, d[2], d[3]) : d[2] > 0 ? (d[1], d[2] - 1, d[3]) : (d[1], d[2], d[3] - 1)
        i = d[1] > 0 ? 1 : d[2] > 0 ? 2 : 3
        deriv(getd(base), i)
    end
    for (der, mono, coef) in PEKERIS_OP
        c = BigFloat(coef(epsv, E, Z, lam))
        iszero(c) && continue
        for (k, x) in getd(der)
            addto!(out, (k[1] + mono[1], k[2] + mono[2], k[3] + mono[3], k[4], k[5]), c * x)
        end
    end
    out
end

struct Moments
    fact::Vector{BigFloat}
    γ::BigFloat
    hs::Matrix{BigFloat}  # hs[k, n] = sum_{j<n} j^-k, k = 1..4
    ζ::Vector{BigFloat}   # ζ(2..4) at index 2..4
end
# Euler–Maclaurin zeta for s >= 2 (accuracy ~1e-60 with M = 1000)
function zeta_em(s::Int; M = 1000)
    B = [1 // 6, -1 // 30, 1 // 42, -1 // 30, 5 // 66, -691 // 2730, 7 // 6]
    z = sum(BigFloat(k)^(-s) for k = 1:M-1) + BigFloat(M)^(1 - s) / (s - 1) + BigFloat(M)^(-s) / 2
    rising = BigFloat(s)
    for (j, b) in enumerate(B)
        z += BigFloat(b) / factorial(big(2j)) * rising * BigFloat(M)^(-s - 2j + 1)
        rising *= (s + 2j - 1) * (s + 2j)
    end
    z
end
function Moments(nmax)
    fact = [factorial(big(n)) for n = 0:nmax]
    hs = zeros(BigFloat, 4, nmax + 1)
    for n = 2:nmax+1, k = 1:4
        hs[k, n] = hs[k, n-1] + 1 / BigFloat(n - 1)^k
    end
    Moments(fact, BigFloat(Base.MathConstants.eulergamma), hs, [zero(BigFloat), zeta_em(2), zeta_em(3), zeta_em(4), zeta_em(5)])
end
function moment(M::Moments, a, b, c, p, q)
    N = a + b + c
    x = N + 3 - q
    @assert x >= 1
    base = M.fact[a+1] * M.fact[b+1] * M.fact[c+1] * M.fact[x] / M.fact[N+3]  # Γ(x) = (x-1)!
    p == 0 && return base
    ψ0 = M.hs[1, x] - M.γ
    p == 1 && return base * ψ0
    ψ1 = M.ζ[2] - M.hs[2, x]
    p == 2 && return base * (ψ0^2 + ψ1)
    ψ2 = -2 * (M.ζ[3] - M.hs[3, x])
    p == 3 && return base * (ψ0^3 + 3ψ0 * ψ1 + ψ2)
    ψ3 = 6 * (M.ζ[4] - M.hs[4, x])
    p == 4 && return base * (ψ0^4 + 6ψ0^2 * ψ1 + 4ψ0 * ψ2 + 3ψ1^2 + ψ3)
    error("p > 4")
end

function inner(M, f::Term, g::Term)
    s = zero(BigFloat)
    for (k1, c1) in f, (k2, c2) in g
        s += c1 * c2 * moment(M, k1[1] + k2[1], k1[2] + k2[2], k1[3] + k2[3], k1[4] + k2[4], k1[5] + k2[5])
    end
    s
end

poly(pairs...) = (t = Term(); for (k, c) in pairs; addto!(t, k, BigFloat(c)); end; t)
symmono(a, b, c, p = 0) = a == b ? poly(((a, b, c, p, 0), 1)) : poly(((a, b, c, p, 0), 1), ((b, a, c, p, 0), 1))
function mul(f::Term, g::Term)
    h = Term()
    for (k1, c1) in f, (k2, c2) in g
        addto!(h, k1 .+ k2, c1 * c2)
    end
    h
end

# r1^2 + r2^2 - r12^2 = (w^2 + 2uw + 2vw - 4uv) / (8 eps^2)
rq(epsv) = poly(((0, 0, 2, 0, 0), 1 / (8epsv^2)), ((1, 0, 1, 0, 0), 2 / (8epsv^2)),
    ((0, 1, 1, 0, 0), 2 / (8epsv^2)), ((1, 1, 0, 0, 0), -4 / (8epsv^2)))
const LOGσ = poly(((0, 0, 0, 1, 0), 1))

function basis(N, variant, epsv; Z = 2)
    κ = (BigFloat(π) - 2) / (3BigFloat(π))
    B = Term[]
    for n = 0:N, a = 0:n, b = a:n-a
        c = n - a - b
        f = symmono(a, b, c)
        if n == 0 && variant in (:fixed, :fixedlog, :fixedlog2)
            f = mergewith(+, f, mul(mul(LOGσ, rq(epsv)), poly(((0, 0, 0, 0, 0), -Z * κ))))
        end
        push!(B, f)
    end
    if variant == :freelog
        push!(B, mul(LOGσ, rq(epsv)))
    end
    if variant in (:freelog, :fixedlog, :fixedlog2)
        for n = 3:N-1, a = 0:n, b = a:n-a  # log sector polynomials of degree 3..N-1
            push!(B, symmono(a, b, n - a - b, 1))
        end
    end
    if variant == :fixedlog2
        for n = 4:N-2, a = 0:n, b = a:n-a  # log^2 sector polynomials of degree 4..N-2
            push!(B, symmono(a, b, n - a - b, 2))
        end
    end
    B
end

function solve(N, variant; epsv = big"1.704", Z = 2, lam = 1, shift = big"-2.9037", bits = 512)
    setprecision(BigFloat, bits) do
        e = BigFloat(epsv)
        B = basis(N, variant, e; Z)
        n = length(B)
        M = Moments(4N + 20)
        PB0 = [apply_P(f, e, 0, Z, lam) for f in B]
        QB = [mergewith(-, apply_P(f, e, 1, Z, lam), PB0[j]) for (j, f) in enumerate(B)]
        H = zeros(BigFloat, n, n); S = zeros(BigFloat, n, n)
        Threads.@threads for i = 1:n
            for j = 1:n
                H[i, j] = -inner(M, B[i], PB0[j])
                S[i, j] = inner(M, B[i], QB[j])
            end
        end
        asym = maximum(abs.(H - H')) / maximum(abs.(H))
        asymS = maximum(abs.(S - S')) / maximum(abs.(S))
        H = (H + H') / 2; S = (S + S') / 2
        F = lu(H - BigFloat(shift) * S)
        x = ones(BigFloat, n); E = BigFloat(0)
        for _ = 1:60
            y = F \ (S * x)
            x = y / norm(y)
            Enew = dot(x, H * x) / dot(x, S * x)
            abs(Enew - E) < BigFloat(10)^(-40) && (E = Enew; break)
            E = Enew
        end
        resid = norm(H * x - E * (S * x)) / norm(H * x)
        posdef = isposdef(Symmetric(S))
        (; E, x, n, asym, asymS, resid, posdef, B)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    E0 = big"-2.903724377034119598311"
    variants = (:none, :fixed, :freelog, :fixedlog)
    open("galerkin_log.out", "w") do io
        println(io, "variant N size E E-E0 asymmetry")
        for N = 2:parse(Int, get(ARGS, 1, "10")), v in variants
            t = @elapsed r = solve(N, v)
            line = @sprintf("%-9s %2d %4d  %.16f  %.4e  %.1e  (%.1fs)", v, N, r.n, Float64(r.E), Float64(r.E - E0), Float64(r.asym), t)
            if v == :freelog
                # ratio of free log coefficient to the constant coefficient vs -Z kappa
                κ = (π - 2) / (3π)
                line *= @sprintf("  c_log/c_0 = %.6f (Fock: %.6f)", Float64(r.x[end-count(f -> any(k[4] == 1 for k in keys(f)), r.B)+1] / r.x[1]), -2κ)
            end
            println(line); println(io, line); flush(io)
        end
    end
end
