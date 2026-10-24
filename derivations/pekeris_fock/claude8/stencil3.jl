# Fast exact Laguerre-stencil generator, v2: dense integer tensors.
# Operator: sum_t c_t u^i v^j w^k eps^a PI^b d_u^du d_v^dv d_w^dw. Terms are grouped by (a,b) (eps, PI powers); in each group the
# rational coefficients are scaled to integers by their LCM. Univariate Laguerre actions x^i d^k L_n = sum_shift p_shift(n) L_{n+shift}
# have integer polynomial coefficients; the stencil is accumulated as dense Int128 tensors in (l,m,n) per shift triple.
const DEG = 16  # max polynomial degree + 1 per index variable
using Base.Checked: checked_add, checked_mul

struct Act1
    shifts::Vector{Int}
    polys::Vector{Vector{Int128}}   # coefficient vectors, index 1 = n^0
    dterm::Bool                      # contains a bare derivative term (must cancel)
    dshifts::Vector{Int}
    dpolys::Vector{Vector{Int128}}
end

pz() = zeros(Int128, DEG)
function pshift(p::Vector{Int128}, j::Int)   # p(n + j)
    r = pz()
    for k = 0:DEG-1
        c = p[k+1]; iszero(c) && continue
        for t = 0:k
            r[t+1] += c * binomial(Int128(k), Int128(t)) * Int128(j)^(k - t)
        end
    end
    r
end
function pmul_n(p::Vector{Int128}, a::Int128, b::Int128)  # p(n) * (a n + b)
    r = pz()
    for k = 0:DEG-1
        c = p[k+1]; iszero(c) && continue
        r[k+1] += b * c
        k + 1 < DEG || (iszero(a) || error("degree overflow")); k + 1 < DEG && (r[k+2] += a * c)
    end
    r
end

# dictionary form during construction: (type, shift) => poly in n (type 0 = L, 1 = bare D)
function xmul(d::Dict{Tuple{Int,Int},Vector{Int128}})
    out = Dict{Tuple{Int,Int},Vector{Int128}}()
    add!(k, p) = (out[k] = get(out, k, pz()) .+ p)
    for ((ty, j), c) in d
        # q = n + j: polynomial in n of c * (affine in q)
        if ty == 0
            add!((0, j + 1), pmul_n(c, Int128(-1), Int128(-(j + 1))))     # -(q+1)
            add!((0, j), pmul_n(c, Int128(2), Int128(2j + 1)))            # 2q+1
            add!((0, j - 1), pmul_n(c, Int128(-1), Int128(-j)))           # -q
        else
            add!((0, j), pmul_n(c, Int128(1), Int128(j)))                  # q L_q
            add!((0, j - 1), pmul_n(c, Int128(-1), Int128(-j)))           # -q L_{q-1}
        end
    end
    filter!(p -> any(!iszero, p.second), out)
end
const CACHE = Dict{Tuple{Int,Int},Dict{Tuple{Int,Int},Vector{Int128}}}()
const CLOCK = ReentrantLock()
function action(a::Int, k::Int)
    lock(CLOCK) do
        get!(CACHE, (a, k)) do
            n = pz(); n[2] = 1
            one_ = pz(); one_[1] = 1
            if k == 0
                d = Dict((0, 0) => copy(one_))
            elseif k == 1
                if a >= 1
                    d = Dict((0, 0) => copy(n), (0, -1) => -n); a -= 1
                else
                    d = Dict((1, 0) => copy(one_))
                end
            elseif k == 2
                if a >= 2
                    xl1 = Dict((0, 0) => copy(n), (0, -1) => -n)
                    t = xmul(xl1)
                    for (kk, c) in xl1; t[kk] = get(t, kk, pz()) .- c; end
                    for (kk, c) in xmul(Dict((0, 0) => copy(n))); t[kk] = get(t, kk, pz()) .- c; end
                    d = filter(p -> any(!iszero, p.second), t); a -= 2
                elseif a == 1
                    d = Dict((0, -1) => -n, (1, 0) => -one_); a -= 1
                else
                    error("bare second derivative")
                end
            else
                error("derivative order > 2")
            end
            for _ = 1:a
                d = xmul(d)
            end
            d
        end
    end
end

struct Term3
    d::NTuple{3,Int}
    m::NTuple{5,Int}
    c::Rational{BigInt}
end

"stencil for one (eps,PI) group with integer coefficients: Dict(shift triple / D-flags => DEG^3 Int128 tensor)"
function group_stencil(terms::Vector{Tuple{NTuple{3,Int},NTuple{3,Int},Int128}})
    acc = Dict{NTuple{3,Tuple{Int,Int}},Array{Int128,3}}()
    d1 = maximum(t[2][1] for t in terms) + 1; d2 = maximum(t[2][2] for t in terms) + 1; d3 = maximum(t[2][3] for t in terms) + 1
    for (d, m, c) in terms
        A = action(m[1], d[1]); B = action(m[2], d[2]); C = action(m[3], d[3])
        for (ka, pa) in A, (kb, pb) in B
            for (kc, pc) in C
                key = (ka, kb, kc)
                T = get!(() -> zeros(Int128, d1, d2, d3), acc, key)
                @inbounds for k3 = 1:d3
                    z3 = pc[k3]; iszero(z3) && continue
                    for k2 = 1:d2
                        z2 = pb[k2]; iszero(z2) && continue
                        f = checked_mul(checked_mul(c, z2), z3)
                        for k1 = 1:d1
                            z1 = pa[k1]; iszero(z1) && continue
                            T[k1, k2, k3] = checked_add(T[k1, k2, k3], checked_mul(f, z1))
                        end
                    end
                end
            end
        end
    end
    filter!(p -> any(!iszero, p.second), acc)
    bad = [k for k in keys(acc) if any(x -> x[1] == 1, k)]
    isempty(bad) || error("non-cancelling bare derivative terms")
    acc
end

# Taylor shift matrix: coefficients of p(x - δ) from p(x)
function shiftmat(δ::Int)
    M = zeros(Int128, DEG, DEG)
    for k = 0:DEG-1, t = 0:k
        M[t+1, k+1] = binomial(Int128(k), Int128(t)) * Int128(-δ)^(k - t)
    end
    M
end
const SM = Dict{Int,Matrix{Int128}}()
smat(δ) = get!(() -> shiftmat(δ), SM, δ)

"block recurrence: offset => Dict((i,j,k,a,b) => Rational coefficient), coefficient polynomial in the equation index"
function block_recurrence(terms::Vector{Term3})
    groups = Dict{NTuple{2,Int},Vector{Term3}}()
    for t in terms
        push!(get!(groups, (t.m[4], t.m[5]), Term3[]), t)
    end
    rec = Dict{NTuple{3,Int},Dict{NTuple{5,Int},Rational{BigInt}}}()
    for ((a, b), ts) in groups
        L = lcm([denominator(t.c) for t in ts])
        it = [(t.d, (t.m[1], t.m[2], t.m[3]), Int128(numerator(t.c) * (L ÷ denominator(t.c)))) for t in ts]
        st = group_stencil(it)
        for (key, T) in st
            δ = (key[1][2], key[2][2], key[3][2]); off = (-δ[1], -δ[2], -δ[3])
            # substitute index -> index - δ in each variable
            S1, S2, S3 = smat(δ[1]), smat(δ[2]), smat(δ[3])
            U = zeros(Int128, size(T))
            for k3 = 1:size(T, 3), k2 = 1:size(T, 2), k1 = 1:size(T, 1)
                v = T[k1, k2, k3]; iszero(v) && continue
                for t1 = 1:k1
                    s1 = S1[t1, k1]; iszero(s1) && continue
                    for t2 = 1:k2
                        s2 = S2[t2, k2]; iszero(s2) && continue
                        for t3 = 1:k3
                            s3 = S3[t3, k3]; iszero(s3) && continue
                            U[t1, t2, t3] = checked_add(U[t1, t2, t3], checked_mul(checked_mul(checked_mul(v, s1), s2), s3))
                        end
                    end
                end
            end
            R = get!(rec, off, Dict{NTuple{5,Int},Rational{BigInt}}())
            for I in CartesianIndices(U)
                v = U[I]; iszero(v) && continue
                k = (I[1] - 1, I[2] - 1, I[3] - 1, a, b)
                R[k] = get(R, k, zero(Rational{BigInt})) + Rational{BigInt}(BigInt(v), L)
            end
        end
    end
    for (k, R) in rec
        filter!(p -> !iszero(p.second), R)
    end
    filter!(p -> !isempty(p.second), rec)
end

function read_terms(path)
    blocks = Dict{String,Vector{Term3}}()
    for line in eachline(path)
        key, ds, ms, num, den = split(line, '|')
        push!(get!(blocks, key, Term3[]), Term3(Tuple(parse.(Int, split(ds))), Tuple(parse.(Int, split(ms))), parse(BigInt, num) // parse(BigInt, den)))
    end
    blocks
end

if abspath(PROGRAM_FILE) == @__FILE__
    # count mode: offsets, monomial coefficients and shell range per block (stencils are regenerated on demand, not stored)
    blocks = read_terms(ARGS[1])
    keysort = sort(collect(keys(blocks)))
    stats = Vector{NTuple{4,Int}}(undef, length(keysort))
    t = @elapsed Threads.@threads :dynamic for i in eachindex(keysort)
        r = block_recurrence(blocks[keysort[i]])
        dN = [sum(o) for o in keys(r)]
        stats[i] = (length(r), sum(length(R) for R in values(r)), minimum(dN), maximum(dN))
        r = nothing
        i % 4 == 0 && GC.gc(false)
    end
    println("blocks $(length(keysort)), total offsets $(sum(s[1] for s in stats)), monomial coefficients $(sum(s[2] for s in stats)), shell range [$(minimum(s[3] for s in stats)), $(maximum(s[4] for s in stats))], largest block $(maximum(s[1] for s in stats)) offsets, time $(round(t; digits = 1)) s")
    open(ARGS[2], "w") do io
        for (k, s_) in zip(keysort, stats); println(io, k, "|", join(s_, " ")); end
    end
end
