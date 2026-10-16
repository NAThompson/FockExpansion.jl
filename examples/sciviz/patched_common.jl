# Shared helpers for the patched-representation figures (patched_*.jl).
#
# Source of the construction: Codex workspace
#     ~/Documents/Codex/2026-09-05/how-much-compu/outputs/boundary_matching/EFFECTIVE_*.txt
# Those files are only read here. Every number in the figures is re-derived in Julia
# and checked against the JSON files written by Codex's budget scripts.
include(joinpath(@__DIR__, "common.jl"))
using Printf, QuadGK

const CODEX = expanduser("~/Documents/Codex/2026-09-05/how-much-compu/outputs")
const BM = joinpath(CODEX, "boundary_matching")
const WITNESS = joinpath(CODEX, "a21_global", "runs", "o10_k8_d80_l1", "result.json")
const NCFILE = joinpath(OUT, "patched.nc")

# Categorical colours (dataviz reference palette, light mode), in fixed order.
const C = (
    blue = colorant"#2a78d6",
    orange = colorant"#eb6834",
    aqua = colorant"#1baf7a",
    yellow = colorant"#eda100",
    magenta = colorant"#e87ba4",
    green = colorant"#008300",
    violet = colorant"#6250d6",
    gray = colorant"#8a8985",
    ink = colorant"#0b0b0b",
    ink2 = colorant"#52514e",
)

# Parameters of the construction (EFFECTIVE_* notes).
const A_SHELL = 1 // 2048        # extraction shell [a/2, a]
const R_STAR = 1 // 16384        # Fock ball radius, original (B = 8192)
const R_STAR_SHARP = 1 // 8192   # after the 2026-10-14 sharpening (B = 4096)
const T_HEAT = 1 // 1024         # heat time t*

# ---------------------------------------------------------------------------
# Minimal JSON reading. The files are small, flat and written by json.dumps; we
# only need a few fields, so read them with regular expressions rather than adding
# a dependency to the shared environment.
jsonstring(txt, key) = match(Regex("\"$key\":\\s*\"([^\"]*)\""), txt).captures[1]
jsonnumber(txt, key) =
    parse(Float64, match(Regex("\"$key\":\\s*([-0-9.eE+]+)"), txt).captures[1])
jsonall(txt, key) =
    [m.captures[1] for m in eachmatch(Regex("\"$key\":\\s*\"?([^\",\\n]*)\"?"), txt)]

"""Convert a parsed JSON array literal (Julia `Expr(:vect)`) to nested vectors,
without evaluating anything."""
tovalue(x::Expr) =
    x.head === :vect ? Any[tovalue(a) for a in x.args] : error("unexpected $(x.head)")
tovalue(x::Union{Number,String}) = x
tovalue(x) = error("unexpected token $x")

function jsonarray(txt, key)
    i = findfirst("\"$key\":", txt)
    j = findnext('[', txt, last(i))
    depth = 0
    k = j
    for (n, ch) in pairs(txt[j:end])
        ch == '[' && (depth += 1)
        ch == ']' && (depth -= 1)
        if depth == 0
            k = j + n - 1
            break
        end
    end
    tovalue(Meta.parse(txt[j:k]))
end

# ---------------------------------------------------------------------------
# Codex witness (the residual-certified seed trial of EFFECTIVE_SEED_CERTIFICATE):
# a 151-term correlated logarithmic Hylleraas function
#     v = exp(-Z s + u/2) * sum_i c_i sum_terms w * s^a t^b u^c (log s)^p,
# s = r1 + r2, t = r1 - r2, u = r12 (outputs/a21_global/variational.py).
# Normalized with v(0) = 1, matching psi00 = 1.
struct Witness
    terms::Vector{NTuple{4,Int}}
    weights::Vector{BigFloat}
    energy::BigFloat
    a21::BigFloat
end

function load_witness(path = WITNESS)
    setprecision(BigFloat, 256) do
        txt = read(path, String)
        basis = jsonarray(txt, "basis")
        coeffs = parse.(BigFloat, jsonarray(txt, "coefficients"))
        terms = NTuple{4,Int}[]
        weights = BigFloat[]
        for (f, c) in zip(basis, coeffs), (key, w) in f
            push!(terms, Tuple(Int.(key)))
            push!(weights, c * parse(BigFloat, w))
        end
        Witness(
            terms,
            weights,
            parse(BigFloat, jsonstring(txt, "energy")),
            parse(BigFloat, jsonstring(txt, "a21")),
        )
    end
end

"""Witness polynomial part minus one, sum_i w_i s^a t^b u^c (log s)^p - 1, in type T."""
function witness_poly(W::Witness, s::T, t::T, u::T) where {T}
    ls = log(s)
    acc = zero(T)
    for ((a, b, c, p), w) in zip(W.terms, W.weights)
        (a, b, c, p) == (0, 0, 0, 0) && continue   # the constant 1 (weight exactly 1)
        acc += T(w) * s^a * t^b * u^c * ls^p
    end
    acc
end

"""Witness at hyperspherical (R, α, θ) minus one, computed so that small differences
from the Fock partial sums survive: v - 1 = expm1(x) + exp(x) * (poly - 1)."""
function witness_m1(W::Witness, R, α, θ; T = BigFloat)
    r1 = T(R) * cos(T(α) / 2)
    r2 = T(R) * sin(T(α) / 2)
    r12 = T(R) * sqrt(1 - sin(T(α)) * cos(T(θ)))
    s, t, u = r1 + r2, r1 - r2, r12
    x = -2s + u / 2
    expm1(x) + exp(x) * witness_poly(W, s, t, u)
end
witness(W::Witness, R, α, θ) = 1 + Float64(witness_m1(W, R, α, θ))

"""L2 norm squared of the witness, ∫ v² d³r₁ d³r₂ = 8π² ∫ r₁² r₂² sinθ v² dr₁ dr₂ dθ,
by tensor Gauss–Legendre in Float64 (adequate for an overall scale factor)."""
function witness_norm2(W::Witness; n = 48, rmax = 16.0, panels = (0.0, 1.0, 3.0, 7.0, 16.0))
    x, w = QuadGK.gauss(n)
    rn, rw = Float64[], Float64[]
    for (a, b) in zip(panels[1:(end-1)], panels[2:end])
        append!(rn, (b - a) / 2 .* x .+ (a + b) / 2)
        append!(rw, (b - a) / 2 .* w)
    end
    θn = π / 2 .* x .+ π / 2
    θw = π / 2 .* w
    tf = Float64.(W.weights)
    acc = zeros(length(rn))
    Threads.@threads for i in eachindex(rn)
        loc = 0.0
        r1 = rn[i]
        for j in eachindex(rn), k in eachindex(θn)
            r2 = rn[j]
            θ = θn[k]
            u = sqrt(r1^2 + r2^2 - 2r1 * r2 * cos(θ))
            s, t = r1 + r2, r1 - r2
            ls = log(s)
            p = 0.0
            for ((a, b, c, q), wt) in zip(W.terms, tf)
                p += wt * s^a * t^b * u^c * ls^q
            end
            v = exp(-2s + u / 2) * p
            loc += rw[i] * rw[j] * θw[k] * r1^2 * r2^2 * sin(θ) * v^2
        end
        acc[i] = loc
    end
    8π^2 * sum(acc)
end

# ---------------------------------------------------------------------------
# Truncated Fock series F_K(R, α, θ) - 1 = sum_{1<=k<=K} R^k sum_p (log R)^p ψ_{k,p}.
const FOCK_TERMS = [
    (1, 0, (a, t) -> psi10(a, t; Z = HE.Z)),
    (2, 0, (a, t) -> psi20(a, t; HE...)),
    (2, 1, (a, t) -> psi21(a, t; Z = HE.Z)),
    (3, 0, (a, t) -> psi30(a, t; HE..., rtol = 1e-12)),
    (3, 1, (a, t) -> psi31(a, t; Z = HE.Z)),
]
fock_angular(α, θ) = [(k, p, f(α, θ)) for (k, p, f) in FOCK_TERMS]
fock_m1(coeffs, R, K) =
    sum(R^k * log(R)^p * c for (k, p, c) in coeffs if k <= K; init = zero(R))

# Smooth radial cutoff: χ = 1 for ρ ≤ r/2, 0 for ρ ≥ r, C^∞ in between.
function chi(ρ, r)
    x = (ρ - r / 2) / (r / 2)
    x <= 0 && return 1.0
    x >= 1 && return 0.0
    f(y) = y > 0 ? exp(-1 / y) : 0.0
    f(1 - x) / (f(1 - x) + f(x))
end

# ---------------------------------------------------------------------------
# NetCDF: every patched_*.jl script writes its own group of output/patched.nc. A
# rerun overwrites the values in place (shapes are deterministic).
function ncgroup(f, name; attrib = Pair{String,Any}[])
    mode = isfile(NCFILE) ? "a" : "c"
    NCDataset(NCFILE, mode) do ds
        if mode == "c"
            ds.attrib["title"] = "Patched global representation of the helium ground state (Codex EFFECTIVE_* construction)"
            ds.attrib["source"] = "FockExpansion.jl examples/sciviz/patched_*.jl"
        end
        g = haskey(ds.group, name) ? ds.group[name] : defGroup(ds, name)
        for (k, v) in attrib
            g.attrib[k] = v
        end
        f(g)
    end
    println("updated ", NCFILE, " group /", name)
end

function putvar!(g, name, data, dims; attrib = Pair{String,Any}[])
    for (d, n) in zip(dims, size(data))
        haskey(g.dim, d) || defDim(g, d, n)
    end
    if haskey(g, name)
        v = g[name]
        v[ntuple(_->Colon(), ndims(data))...] = data
        for (k, val) in attrib
            v.attrib[k] = val
        end
    else
        defVar(g, name, data, dims; attrib = attrib)
    end
end
