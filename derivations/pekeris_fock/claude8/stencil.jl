# Fast exact Laguerre-stencil generator (replaces the sympy operator_to_stencil / stencil_to_recurrence of sym/recurrence.py).
# Input: polynomial-coefficient differential operators  sum c * u^i v^j w^k eps^a PI^b * d_u^du d_v^dv d_w^dw  (c rational),
# acting on L_l(u) L_m(v) L_n(w). Output: recurrence  sum_off coeff(l,m,n) A(l+off) = 0, coefficient polynomials in
# (l, m, n, eps, PI) with exact Rational{BigInt} coefficients. Laguerre rules:
#   x L_n = -(n+1) L_{n+1} + (2n+1) L_n - n L_{n-1},  x L_n' = n L_n - n L_{n-1},  x L_n'' = (x-1) L_n' - n L_n.
# A bare L_n' (no x factor) is kept as a 'D' term and must cancel in the total (else the stencil would be infinite).

const Q = Rational{BigInt}
const UPoly = Dict{Int,Q}                   # univariate polynomial in the Laguerre index (power => coeff)
const MKey = NTuple{5,Int}                   # exponents of (l, m, n, eps, PI)
const MPoly = Dict{MKey,Q}

upoly(pairs...) = (p = UPoly(); for (k, c) in pairs; p[k] = get(p, k, zero(Q)) + Q(c); end; p)
function umul(a::UPoly, b::UPoly)
    r = UPoly()
    for (i, x) in a, (j, y) in b
        r[i+j] = get(r, i + j, zero(Q)) + x * y
    end
    filter!(p -> !iszero(p.second), r)
end
uadd!(a::UPoly, b::UPoly, s = one(Q)) = (for (k, x) in b; a[k] = get(a, k, zero(Q)) + s * x; end; filter!(p -> !iszero(p.second), a))
# polynomial in n shifted: p(n + j)
function ushift(p::UPoly, j::Int)
    r = UPoly()
    for (k, c) in p, t = 0:k
        r[t] = get(r, t, zero(Q)) + c * binomial(BigInt(k), BigInt(t)) * BigInt(j)^(k - t)
    end
    filter!(q -> !iszero(q.second), r)
end

# 1D action dictionaries: (type, shift) => polynomial in n, type 0 = L, 1 = D (bare derivative)
const Act = Dict{Tuple{Int,Int},UPoly}
function xmul(d::Act)
    out = Act()
    add!(k, p) = (out[k] = uadd!(get(out, k, UPoly()), p))
    for ((ty, j), c) in d
        q = ushift(upoly(1 => 1), j)                        # q = n + j
        if ty == 0
            add!((0, j + 1), umul(c, uadd!(ushift(upoly(1 => -1), j), upoly(0 => -1))))   # -(q+1)
            add!((0, j), umul(c, uadd!(ushift(upoly(1 => 2), j), upoly(0 => 1))))         # 2q+1
            add!((0, j - 1), umul(c, ushift(upoly(1 => -1), j)))                            # -q
        else
            add!((0, j), umul(c, q)); add!((0, j - 1), umul(c, ushift(upoly(1 => -1), j)))
        end
    end
    filter!(p -> !isempty(p.second), out)
end
const ACT_CACHE = Dict{Tuple{Int,Int},Act}()
function action1d(a::Int, k::Int)
    get!(ACT_CACHE, (a, k)) do
        n = upoly(1 => 1)
        if k == 0
            d = Act((0, 0) => upoly(0 => 1))
        elseif k == 1
            if a >= 1
                d = Act((0, 0) => n, (0, -1) => upoly(1 => -1)); a -= 1
            else
                d = Act((1, 0) => upoly(0 => 1))
            end
        elseif k == 2
            if a >= 2
                xl1 = Act((0, 0) => n, (0, -1) => upoly(1 => -1))      # x L'
                t = xmul(xl1)
                for (kk, c) in xl1; t[kk] = uadd!(get(t, kk, UPoly()), c, -one(Q)); end
                t2 = xmul(Act((0, 0) => n))
                for (kk, c) in t2; t[kk] = uadd!(get(t, kk, UPoly()), c, -one(Q)); end
                d = filter(p -> !isempty(p.second), t); a -= 2
            elseif a == 1
                d = Act((0, -1) => upoly(1 => -1), (1, 0) => upoly(0 => -1)); a -= 1
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

# operator terms: (deriv, mono(i,j,k,a,b), coeff)
struct Term3
    d::NTuple{3,Int}
    m::NTuple{5,Int}
    c::Q
end

"Stencil keyed by the three (type, shift) pairs, coefficient polynomial in (l,m,n,eps,PI) evaluated at the SOURCE index."
function stencil(terms::Vector{Term3})
    st = Dict{NTuple{3,Tuple{Int,Int}},MPoly}()
    for t in terms
        A = action1d(t.m[1], t.d[1]); B = action1d(t.m[2], t.d[2]); C = action1d(t.m[3], t.d[3])
        for (ka, pa) in A, (kb, pb) in B, (kc, pc) in C
            key = (ka, kb, kc)
            acc = get!(st, key, MPoly())
            for (ia, ca) in pa, (ib, cb) in pb, (ic, cc) in pc
                mk = (ia, ib, ic, t.m[4], t.m[5])
                acc[mk] = get(acc, mk, zero(Q)) + t.c * ca * cb * cc
            end
        end
    end
    for (k, p) in st
        filter!(q -> !iszero(q.second), p)
    end
    filter!(p -> !isempty(p.second), st)
    bad = [k for k in keys(st) if any(x -> x[1] == 1, k)]
    isempty(bad) || error("non-cancelling bare derivative terms: $(bad[1:min(3,end)])")
    st
end

"recurrence: offset (= -shift) => coefficient polynomial in the EQUATION index (l -> l - shift)"
function recurrence(st)
    rec = Dict{NTuple{3,Int},MPoly}()
    for (key, p) in st
        δ = (key[1][2], key[2][2], key[3][2])
        off = (-δ[1], -δ[2], -δ[3])
        acc = get!(rec, off, MPoly())
        for ((i, j, k, a, b), c) in p
            # (l - δ1)^i (m - δ2)^j (n - δ3)^k
            for ti = 0:i, tj = 0:j, tk = 0:k
                cf = c * binomial(BigInt(i), BigInt(ti)) * BigInt(-δ[1])^(i - ti) *
                     binomial(BigInt(j), BigInt(tj)) * BigInt(-δ[2])^(j - tj) *
                     binomial(BigInt(k), BigInt(tk)) * BigInt(-δ[3])^(k - tk)
                mk = (ti, tj, tk, a, b)
                acc[mk] = get(acc, mk, zero(Q)) + cf
            end
        end
    end
    for (k, p) in rec
        filter!(q -> !iszero(q.second), p)
    end
    filter!(p -> !isempty(p.second), rec)
end

function read_terms(path)
    blocks = Dict{String,Vector{Term3}}()
    for line in eachline(path)
        key, ds, ms, num, den = split(line, '|')
        d = Tuple(parse.(Int, split(ds))); m = Tuple(parse.(Int, split(ms)))
        push!(get!(blocks, key, Term3[]), Term3(d, m, Q(parse(BigInt, num), parse(BigInt, den))))
    end
    blocks
end

function write_rec(io, key, rec)
    for (off, p) in rec, (mk, c) in p
        println(io, key, "|", join(off, " "), "|", join(mk, " "), "|", numerator(c), "|", denominator(c))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    blocks = read_terms(ARGS[1])
    keysort = sort(collect(keys(blocks)))
    results = Vector{Any}(undef, length(keysort))
    t = @elapsed Threads.@threads for i in eachindex(keysort)
        results[i] = recurrence(stencil(blocks[keysort[i]]))
    end
    total = sum(length(r) for r in results)
    println("blocks $(length(keysort)), total offsets $total, time $(round(t; digits = 1)) s")
    open(ARGS[2], "w") do io
        for (k, r) in zip(keysort, results)
            write_rec(io, k, r)
        end
    end
end
