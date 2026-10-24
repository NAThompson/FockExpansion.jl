# Phase 2 C: Rayleigh–Ritz with Fock "atoms" in (s, x, y) coordinates,
#   s = r1 + r2, x = (r1 - r2)/s, y = r12/s,  Ψ = exp(-ε s) Σ c s^n log(s)^p x^a y^b g(x, y).
# Radial integrals are exact (Γ and polygamma); angular integrals use a tanh–sinh rule on the half
# triangle 0 ≤ x ≤ y ≤ 1 (Duffy map x = y t), in BigFloat. Entries are assembled from moment tables
#   M[h, h'](a, b) = ∫ x^a y^b h h' dx dy,  h ∈ {g, ∂x g, ∂y g} for every atom g.
using LinearAlgebra, Printf, ForwardDiff
using FockExpansion
import ClausenFunctions

# ForwardDiff rule for Cl₂ (d/dt Cl₂(t) = -log|2 sin(t/2)|), so the repo's chi can be differentiated.
function ClausenFunctions.cl2(d::ForwardDiff.Dual{T}) where {T}
    t = ForwardDiff.value(d)
    ForwardDiff.Dual{T}(ClausenFunctions.cl2(t), -log(abs(2sin(t / 2))) * ForwardDiff.partials(d))
end

const Z = 2
const E0 = big"-2.903724377034119598311159245194404446696925310"
const A21 = big"0.47674787900"

# ---------------- angular quadrature ----------------
function tanhsinh01(h, kmax)
    pts = BigFloat[]; wts = BigFloat[]; comp = BigFloat[]  # comp = 1 - t
    for k = -kmax:kmax
        τ = k * h
        q = BigFloat(π) / 2 * sinh(τ)
        e = exp(-2abs(q))
        # t = (1 + tanh q)/2 ; 1 - t = 1/(1 + exp(2q))
        if q >= 0
            onem = e / (1 + e); t = 1 - onem
        else
            t = e / (1 + e); onem = 1 - t
        end
        w = h * BigFloat(π) / 2 * cosh(τ) / (2cosh(q)^2)
        push!(pts, t); push!(wts, w); push!(comp, onem)
    end
    pts, wts, comp
end

struct Nodes
    x::Vector{BigFloat}
    y::Vector{BigFloat}
    w::Vector{BigFloat}
end
function nodes(h, kmax)
    t, wt, _ = tanhsinh01(h, kmax)
    X = BigFloat[]; Y = BigFloat[]; W = BigFloat[]
    cut = BigFloat(10)^(-28)  # drop nodes whose total weight is negligible but whose μ = ±1 would round over
    for i in eachindex(t), j in eachindex(t)
        y = t[j]; x = y * t[i]
        (y < cut || 1 - t[i] < cut || 1 - y < cut) && continue
        # 1 ∓ μ = 2(y² - x²)/(1 - x²) and 2(1 - y²)/(1 - x²): keep them resolvable at this precision
        (2 * (y - x) * (y + x) / (1 - x^2) < BigFloat(10)^(-50) || 2 * (1 - y^2) / (1 - x^2) < BigFloat(10)^(-50)) && continue
        push!(X, x); push!(Y, y); push!(W, wt[i] * wt[j] * y)
    end
    Nodes(X, Y, W)
end

# ---------------- atoms as functions of (x, y) at s = 1 ----------------
function angles(x, y)
    r1 = (1 + x) / 2
    r2 = (1 - x) / 2
    2atan(r2 / r1), (r1^2 + r2^2 - y^2) / (2r1 * r2)
end
# θ without cancellation: 1 - cos θ = 2(y² - x²)/(1 - x²)
theta(x, y) = 2asin(sqrt((y - x) * (y + x) / (1 - x^2)))
function chi_xy(x, y)
    α, μ = angles(x, y)
    FockExpansion.chi(α, μ)
end
lR_xy(x, y) = log((1 + x^2) / 2) / 2

# ---- Codex's minimal order-2 atoms (global letter forms, s = 1; d = x, r12 = y) ----
function o2parts(x, y)
    r1 = (1 + x) / 2; r2 = (1 - x) / 2; R2 = (1 + x^2) / 2
    P = sqrt(2R2 - y^2)
    α = 2atan(r2 / r1); β = asin((R2 - y^2) / R2)
    L5 = y * P + x
    p2 = 2r1 * r2
    cl = ClausenFunctions.cl2
    DL = cl(α - β) + cl(oftype(x, π) + β - α) - cl(oftype(x, π) - α - β) - cl(α + β)
    Ω = β * log(p2^2 / L5^2) + (β + α) * log((1 - y^2) / p2) + (β - α) * log((y - x) * (y + x) / p2) - DL
    lb = log((x + y)^2 / L5)
    (; r1, r2, R2, P, α, β, Ω, lb)
end
U_xy(x, y) = (q = o2parts(x, y); q.Ω + oftype(x, π) * q.lb)                 # exchange-odd
V_xy(x, y) = (q = o2parts(x, y); q.P / y * (q.β + oftype(x, π) / 2))       # exchange-even
alp_xy(x, y) = 2atan((1 - x) / (1 + x)) - oftype(x, π) / 2                 # exchange-odd
Lt_xy(x, y) = log((1 + y) / sqrt((1 + x^2) / 2)) + (oftype(x, π) - 2) / oftype(x, π) * log(sqrt((1 + x^2) / 2))  # L - ((pi-2)/pi) log(s/R)
# chi from Codex's formula (1), for a consistency check of U, V, L against the repo's chi
function chi_codex(x, y)
    q = o2parts(x, y); R2 = q.R2
    s = 2q.r1 * q.r2 / R2; c = x / R2; v = (R2 - y^2) / R2; g = y * q.P / R2
    L = log((1 + y) / sqrt(R2)); U = q.Ω + oftype(x, π) * q.lb
    h0 = -(3oftype(x, π) + 10 - 16 * oftype(x, Base.MathConstants.catalan)) / (24oftype(x, π))
    (g + s - 2y / R2) / 6 - v * L / 3 + (c * U + 2g * q.β) / (6oftype(x, π)) + h0 * v
end
function psi30h_xy(x, y)
    α, μ = angles(x, y)
    hx = sqrt((1 + x^2) / 2)
    hx^3 * psi30(α, theta(x, y); Z = 2.0, E = Float64(E0), a21 = Float64(A21), rtol = 1e-13)
end

const PSI30_DIGITS = Ref(35)
function psi30h_big(x, y)
    α, μ = angles(x, y)
    hx = sqrt((1 + x^2) / 2)
    hx^3 * psi30(α, theta(x, y); Z = big(2.0), E = BigFloat(E0), a21 = BigFloat(A21), rtol = BigFloat(10)^(-PSI30_DIGITS[]))
end

"value and gradient of an atom at all nodes (BigFloat via ForwardDiff, or Float64 for psi30)"
function atom_table(f, nd::Nodes; float64 = false)
    n = length(nd.x)
    val = zeros(BigFloat, n); gx = zeros(BigFloat, n); gy = zeros(BigFloat, n)
    Threads.@threads for q = 1:n
        g = ForwardDiff.gradient(z -> f(z[1], z[2]), [nd.x[q], nd.y[q]])
        val[q] = f(nd.x[q], nd.y[q]); gx[q] = g[1]; gy[q] = g[2]
    end
    if false
        nbad = count(q -> !isfinite(gx[q]) || !isfinite(gy[q]), 1:n)
        for q = 1:n
            isfinite(gx[q]) || (gx[q] = 0); isfinite(gy[q]) || (gy[q] = 0)
        end
        @info "psi30 atom: zeroed non-finite gradients at $nbad nodes (max weight $(Float64(maximum(nd.w[q] for q in 1:n if gx[q] == 0 && gy[q] == 0; init = big(0)))))"
    end
    val, gx, gy
end

# ---------------- moment tables ----------------
# h-functions are identified by (atom, deriv) with deriv ∈ (:v, :x, :y); the constant atom :one has only :v.
struct Moments
    tab::Dict{Tuple{Symbol,Symbol,Symbol,Symbol},Matrix{BigFloat}}  # (g1, d1, g2, d2) -> M[a+1, b+1]
end
function moments(nd::Nodes, atoms::Dict{Symbol,NTuple{3,Vector{BigFloat}}}, amax, bmax)
    n = length(nd.x)
    Xp = [nd.x[q]^a for q = 1:n, a = 0:amax]
    Yp = [nd.y[q]^b for q = 1:n, b = 0:bmax]
    hs = Tuple{Symbol,Symbol,Vector{BigFloat}}[]
    for (g, (v, gx, gy)) in atoms
        push!(hs, (g, :v, v))
        g == :one || (push!(hs, (g, :x, gx)); push!(hs, (g, :y, gy)))
    end
    tab = Dict{Tuple{Symbol,Symbol,Symbol,Symbol},Matrix{BigFloat}}()
    lk = ReentrantLock()
    pairs = [(i, j) for i in eachindex(hs) for j in i:length(hs)]
    Threads.@threads for (i, j) in pairs
        g1, d1, v1 = hs[i]; g2, d2, v2 = hs[j]
        wv = nd.w .* v1 .* v2 .* 2  # factor 2: half triangle, symmetric total integrands
        M = transpose(Xp) * (wv .* Yp)
        lock(lk) do
            tab[(g1, d1, g2, d2)] = M
            tab[(g2, d2, g1, d1)] = M
        end
    end
    Moments(tab)
end

# ---------------- basis functions ----------------
# A term: s^n log^p s × Σ_k c_k x^a y^b g  (angular part as Dict{(a,b,g)=>c}); a basis function is a vector of terms.
const Ang = Dict{Tuple{Int,Int,Symbol},BigFloat}
struct Term
    n::Int
    p::Int
    Y::Ang
end
addang!(d::Ang, k, c) = (d[k] = get(d, k, zero(BigFloat)) + c; d)
# angular derivative of x^a y^b g: returns Dict over (a, b, g, deriv)
function dang(Y::Ang, var)
    out = Dict{Tuple{Int,Int,Symbol,Symbol},BigFloat}()
    for ((a, b, g), c) in Y
        k = var == :x ? a : b
        if k > 0
            key = var == :x ? (a - 1, b, g, :v) : (a, b - 1, g, :v)
            out[key] = get(out, key, zero(BigFloat)) + c * k
        end
        if g != :one
            key = (a, b, g, var)
            out[key] = get(out, key, zero(BigFloat)) + c
        end
    end
    out
end
plain(Y::Ang) = Dict((a, b, g, :v) => c for ((a, b, g), c) in Y)
# multiply angular dict by polynomial Σ c x^a y^b (Dict{(a,b)=>c})
function polymul(d, poly)
    out = Dict{Tuple{Int,Int,Symbol,Symbol},BigFloat}()
    for ((a, b, g, dv), c) in d, ((pa, pb), pc) in poly
        key = (a + pa, b + pb, g, dv)
        out[key] = get(out, key, zero(BigFloat)) + c * pc
    end
    out
end
function angint(M::Moments, d1, d2)
    s = zero(BigFloat)
    for ((a1, b1, g1, v1), c1) in d1, ((a2, b2, g2, v2), c2) in d2
        s += c1 * c2 * M.tab[(g1, v1, g2, v2)][a1+a2+1, b1+b2+1]
    end
    s
end

# radial ∫ e^{-2εs} s^m log^q s ds = Γ(m+1)(2ε)^{-(m+1)} × Bell_q(ψ0 - log 2ε, ψ1, ψ2, ψ3)
struct Radial
    ε::BigFloat
    cache::Dict{Tuple{Int,Int},BigFloat}
end
Radial(ε) = Radial(ε, Dict{Tuple{Int,Int},BigFloat}())
function polygamma_int(k, x::Int)  # ψ^{(k)}(x) for integer x ≥ 1
    γ = BigFloat(Base.MathConstants.eulergamma)
    k == 0 && return sum((1 / BigFloat(j) for j = 1:x-1); init = zero(BigFloat)) - γ
    ζ = k == 1 ? BigFloat(π)^2 / 6 : k == 3 ? BigFloat(π)^4 / 90 : zeta3()
    (-1)^(k + 1) * factorial(k) * (ζ - sum((1 / BigFloat(j)^(k + 1) for j = 1:x-1); init = zero(BigFloat)))
end
# ζ(3) = (5/2) Σ (-1)^(k+1) / (k³ binomial(2k, k)), geometric convergence
zeta3() = 5 * sum((-1)^(k + 1) / (BigFloat(k)^3 * binomial(big(2k), big(k))) for k = 1:200) / 2
function rad(R::Radial, m, q)
    get!(R.cache, (m, q)) do
        a = 2R.ε
        base = factorial(big(m)) / a^(m + 1)
        q == 0 && return base
        p0 = polygamma_int(0, m + 1) - log(a)
        q == 1 && return base * p0
        p1 = polygamma_int(1, m + 1)
        q == 2 && return base * (p0^2 + p1)
        p2 = polygamma_int(2, m + 1)
        q == 3 && return base * (p0^3 + 3p0 * p1 + p2)
        p3 = polygamma_int(3, m + 1)
        q == 4 && return base * (p0^4 + 6p0^2 * p1 + 4p0 * p2 + 3p1^2 + p3)
        error("q > 4")
    end
end

# derivative pieces of one term: list of (m-shift, q, coefficient multiplier, angular dict) for ∂r1, ∂r2, ∂r12
function deriv_terms(t::Term, ε)
    Yx = dang(t.Y, :x); Yy = dang(t.Y, :y); Y = plain(t.Y)
    one_m_x = Dict((0, 0) => big(1), (1, 0) => big(-1)); one_p_x = Dict((0, 0) => big(1), (1, 0) => big(1)); yv = Dict((0, 1) => big(1))
    function comb(parts...)
        out = Dict{Tuple{Int,Int,Symbol,Symbol},BigFloat}()
        for (d, c) in parts, (k, v) in d
            out[k] = get(out, k, zero(BigFloat)) + c * v
        end
        out
    end
    nY = comb((Y, BigFloat(t.n)))
    A1 = comb((nY, 1), (polymul(Yx, one_m_x), 1), (polymul(Yy, yv), -1))
    A2 = comb((nY, 1), (polymul(Yx, one_p_x), -1), (polymul(Yy, yv), -1))
    d1 = [(-1, t.p, A1), (0, t.p, comb((Y, -ε)))]
    d2 = [(-1, t.p, A2), (0, t.p, comb((Y, -ε)))]
    if t.p > 0
        push!(d1, (-1, t.p - 1, comb((Y, BigFloat(t.p))))); push!(d2, (-1, t.p - 1, comb((Y, BigFloat(t.p)))))
    end
    d3 = [(-1, t.p, Yy)]
    d1, d2, d3
end

# angular weights (after multiplying by the measure (1-x²) y / 8)
const W0 = Dict((0, 1) => big(1) / 8, (2, 1) => -big(1) / 8)                       # measure
const W11 = W0; const W33 = Dict(k => 2v for (k, v) in W0)
const W13 = Dict((1, 0) => big(1) / 8, (0, 2) => big(1) / 8, (2, 0) => -big(1) / 8, (1, 2) => -big(1) / 8)  # (A/2)·w0 = (1-x)(x+y²)/8
const W23 = Dict((1, 0) => -big(1) / 8, (0, 2) => big(1) / 8, (2, 0) => -big(1) / 8, (1, 2) => big(1) / 8)  # (B/2)·w0 = (1+x)(y²-x)/8
const WV = Dict((0, 1) => -big(Z) / 2, (0, 0) => big(1) / 8, (2, 0) => -big(1) / 8)  # s·V·w0

function matrices(basis::Vector{Vector{Term}}, M::Moments, R::Radial, ε)
    nb = length(basis)
    H = zeros(BigFloat, nb, nb); S = zeros(BigFloat, nb, nb)
    D = [[deriv_terms(t, ε) for t in f] for f in basis]
    nmax = maximum(t.n for f in basis for t in f)
    for m = 0:2nmax+8, q = 0:4  # fill the radial cache before threading (Dict is not thread-safe)
        rad(R, m, q)
    end
    Threads.@threads for i = 1:nb
        for j = i:nb
            h = zero(BigFloat); s = zero(BigFloat)
            for (ti, t1) in enumerate(basis[i]), (tj, t2) in enumerate(basis[j])
                p1, p2 = plain(t1.Y), plain(t2.Y)
                s += rad(R, t1.n + t2.n + 5, t1.p + t2.p) * angint(M, p1, polymul(p2, W0))
                h += rad(R, t1.n + t2.n + 4, t1.p + t2.p) * angint(M, p1, polymul(p2, WV))
                Da, Db = D[i][ti], D[j][tj]
                for (wk, (A, B)) in ((W11, (1, 1)), (W11, (2, 2)), (W33, (3, 3)), (W13, (1, 3)), (W13, (3, 1)), (W23, (2, 3)), (W23, (3, 2)))
                    for (m1, q1, a1) in Da[A], (m2, q2, a2) in Db[B]
                        h += rad(R, t1.n + t2.n + 5 + m1 + m2, q1 + q2) * angint(M, a1, polymul(a2, wk)) / 2
                    end
                end
            end
            H[i, j] = h; H[j, i] = h; S[i, j] = s; S[j, i] = s
        end
    end
    H, S
end

function ground(H, S; shift = big"-2.9037244")
    F = lu(H - shift * S)
    x = ones(BigFloat, size(H, 1)); E = zero(BigFloat)
    for _ = 1:80
        y = F \ (S * x); x = y / norm(y)
        En = dot(x, H * x) / dot(x, S * x)
        abs(En - E) < BigFloat(10)^(-45) && (E = En; break)
        E = En
    end
    E, x
end

# Gram condition number estimate (power iteration on S and on S⁻¹ via Cholesky)
function condest(S)
    C = cholesky(Symmetric(S))
    x = ones(BigFloat, size(S, 1)); λmax = zero(BigFloat); λmin = zero(BigFloat)
    for _ = 1:200
        y = S * x; λmax = norm(y) / norm(x); x = y / norm(y)
    end
    x = ones(BigFloat, size(S, 1))
    for _ = 1:200
        y = C \ x; λmin = norm(x) / norm(y); x = y / norm(y)
    end
    λmax / λmin
end

# ---------------- basis ----------------
const κ = (BigFloat(π) - 2) / (3BigFloat(π))
rq() = Ang((0, 0, :one) => big(1) / 2, (2, 0, :one) => big(1) / 2, (0, 2, :one) => -big(1))  # (r1²+r2²-r12²)/s²
scaled(Y::Ang, c) = Ang(k => c * v for (k, v) in Y)
mono(a, b, g = :one) = Ang((a, b, g) => big(1))
function build(N, variant; ε, Na = N - 2)
    B = Vector{Term}[]
    fock0 = [Term(0, 0, mono(0, 0)), Term(2, 1, scaled(rq(), -Z * κ)), Term(2, 0, scaled(rq(), -Z * κ * log(2ε)))]
    push!(B, variant == :none ? [Term(0, 0, mono(0, 0))] : fock0)
    for n = 1:N, a = 0:2:n, b = 0:n-a
        push!(B, [Term(n, 0, mono(a, b))])
    end
    variant == :none && return B
    for n = 3:N-1, a = 0:2:n, b = 0:n-a
        push!(B, [Term(n, 1, mono(a, b))])
    end
    variant == :log1 && return B
    for n = 4:N-2, a = 0:2:n, b = 0:n-a
        push!(B, [Term(n, 2, mono(a, b))])
    end
    variant == :log2 && return B
    for g in (:chi, :lR), n = 2:Na, a = 0:2:n, b = 0:n-a
        push!(B, [Term(n, 0, mono(a, b, g))])
    end
    variant == :atoms2 && return B
    for n = 3:Na, a = 0:2:n-3, b = 0:n-3-a
        push!(B, [Term(n, 0, mono(a, b, :psi30h))])
    end
    for n = 4:Na, a = 0:2:n-4, b = 0:n-4-a
        push!(B, [Term(n, 1, mono(a, b, :lR))]); push!(B, [Term(n, 0, mono(a, b, :lR2))])
    end
    variant == :atoms4 && return B
    error("unknown variant")
end
const ODD = Set([:U, :alp])
function build3(N, variant; ε, Na = N - 2)
    B = build(N, :log2; ε)
    atoms = variant in (:min5, :min5psi) ? (:Lt, :alp, :V, :U) : (:chi, :lR)
    for g in atoms, n = 2:Na, a = (g in ODD ? 1 : 0):2:n, b = (g == :V ? 1 : 0):n-a  # V ~ pi P/r12 at r12 -> 0: require a factor y (Codex pole audit)
        push!(B, [Term(n, 0, mono(a, b, g))])
    end
    if variant in (:min5psi, :atoms2psi)
        for n = 3:Na, a = 0:2:n-3, b = 0:n-3-a
            push!(B, [Term(n, 0, mono(a, b, :psi30h))])
        end
    end
    B
end
