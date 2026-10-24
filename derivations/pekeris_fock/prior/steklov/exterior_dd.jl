# High-precision exterior Ritz assembly.
#
# The Hylleraas correction block is integrated in r12 analytically. Every integrand is a
# polynomial in r12 times a function of (r1, r2), because the S-state measure
# r1 r2 r12 cancels 1/r12 and the direction cosines. A 2D Gauss rule in (R, phi) is left over.
# Lift blocks contain exp(g r12/2) and use a 3D Gauss rule in (R, phi, u).
# All arithmetic is in the element type T (Double64 by default). Without it, the Gram-type
# matrices lose the near-dependent directions of the Hylleraas basis
# (see ext_test2.jl: Float64 results move by 1e-5 with the truncation level).

using DoubleFloats

# Radial rule on [Rs, inf): Gauss–Legendre panels with breakpoints Rs + edges, then a
# Gauss–Laguerre tail. Weights include exp(t) so plain integrands can be summed.
# edges = nothing gives the pure Laguerre rule R = Rs + t/beta.
function radial_rule(Rs, ::Type{T}; nR = 60, beta = 2.5, edges = nothing, npanel = 16) where {T}
    if edges === nothing
        t, wt = gauss_laguerre(nR, T)
        return T(Rs) .+ t ./ T(beta), wt .* exp.(t) ./ T(beta)
    end
    Rn = T[]
    wn = T[]
    a = T(Rs)
    for e in edges
        b = T(Rs) + T(e)
        x, w = gauss_legendre(npanel, a, b, T)
        append!(Rn, x)
        append!(wn, w)
        a = b
    end
    t, wt = gauss_laguerre(nR, T)
    append!(Rn, a .+ t ./ T(beta))
    append!(wn, wt .* exp.(t) ./ T(beta))
    return Rn, wn
end

# 2D rule: nodes (r1, r2) and weights for dr1 dr2 with the (8/pi) r1 r2 factor of the
# measure, doubled for exchange (phi in [0, pi/4]). The r12 factor goes into I(p).
function exterior_rule_2d(Rs, ::Type{T}; nR = 60, nphi = 40, beta = 2.5, edges = nothing,
                          npanel = 16) where {T}
    Rn, wn = radial_rule(Rs, T; nR, beta, edges, npanel)
    ph, wph = gauss_legendre(nphi, zero(T), T(pi) / 4, T)
    r1 = T[]
    r2 = T[]
    w = T[]
    for (R, wR) in zip(Rn, wn), (p, wp) in zip(ph, wph)
        c, s = cos(p), sin(p)
        push!(r1, R * c)
        push!(r2, R * s)
        push!(w, 2 * T(8) / T(pi) * (R * c) * (R * s) * R * wR * wp)
    end
    return (; r1, r2, w)
end

function exterior_rule_3d(Rs, ::Type{T}; nR = 60, nphi = 40, nu = 30, beta = 2.5,
                          edges = nothing, npanel = 16) where {T}
    Rn, wn = radial_rule(Rs, T; nR, beta, edges, npanel)
    ph, wph = gauss_legendre(nphi, zero(T), T(pi) / 4, T)
    xu, wu = gauss_legendre(nu, T)
    r1 = T[]
    r2 = T[]
    r12 = T[]
    w = T[]
    for (R, wR0) in zip(Rn, wn), (p, wp) in zip(ph, wph)
        wR = wR0 * R^5
        c, s = cos(p), sin(p)
        lo, hi = c - s, c + s
        for (x, wx) in zip(xu, wu)
            uu = (hi - lo) / 2 * x + (hi + lo) / 2
            push!(r1, R * c)
            push!(r2, R * s)
            push!(r12, R * uu)
            push!(w, 2 * T(8) / T(pi) * c * s * uu * wp * wx * (hi - lo) / 2 * wR)
        end
    end
    return (; r1, r2, r12, w)
end

# phi_b(r1, r2) = (R - Rs) [p(r1, r2) + p(r2, r1)] and its r1, r2 derivatives,
# with p = r1^i r2^j exp(-za r1 - zb r2). The full term is phi_b * r12^k.
function radial_part(t::Term, r1::T, r2::T, Rs) where {T}
    R = sqrt(r1^2 + r2^2)
    h = R - T(Rs)
    za, zb = T(t.za), T(t.zb)
    p = r1^t.i * r2^t.j * exp(-za * r1 - zb * r2)
    q = r2^t.i * r1^t.j * exp(-za * r2 - zb * r1)
    p1 = (t.i / r1 - za) * p
    p2 = (t.j / r2 - zb) * p
    q1 = (t.j / r1 - zb) * q
    q2 = (t.i / r2 - za) * q
    v = h * (p + q)
    d1 = r1 / R * (p + q) + h * (p1 + q1)
    d2 = r2 / R * (p + q) + h * (p2 + q2)
    return v, d1, d2
end

function correction_block(terms, rule2, Rs, E, Z, g, ::Type{T}) where {T}
    nb = length(terms)
    all(t -> t.zc == 0, terms) || throw(ArgumentError("analytic r12 needs zc = 0"))
    kvec = [t.k for t in terms]
    pmax = 2maximum(kvec) + 2
    nn = length(rule2.w)
    nt = Threads.nthreads()
    partial = [zeros(T, nb, nb) for _ = 1:nt]
    Threads.@threads for chunk = 1:nt
        Acc = partial[chunk]
        v = zeros(T, nb)
        d1 = zeros(T, nb)
        d2 = zeros(T, nb)
        I = zeros(T, pmax + 2)        # I[p + 2] = int r12^p dr12, p >= -1
        for node = chunk:nt:nn
            x1, x2, w = rule2.r1[node], rule2.r2[node], rule2.w[node]
            lo, hi = abs(x1 - x2), x1 + x2
            I[1] = zero(T)            # p = -1 only ever multiplies zero
            for p = 0:pmax
                I[p+2] = (hi^(p + 1) - lo^(p + 1)) / (p + 1)
            end
            for b = 1:nb
                v[b], d1[b], d2[b] = radial_part(terms[b], x1, x2, Rs)
            end
            Vn = 2 * (-T(Z) / x1 - T(Z) / x2 - T(E))
            dd = x1^2 - x2^2
            inv1 = 1 / (2x1)
            inv2 = 1 / (2x2)
            for b = 1:nb, c = b:nb
                k, kp = kvec[b], kvec[c]
                Ip1 = I[k+kp+3]
                Ip0 = I[k+kp+2]
                Im1 = I[k+kp+1]
                vv = v[b] * v[c]
                val = (d1[b] * d1[c] + d2[b] * d2[c] + Vn * vv) * Ip1 + 2 * T(g) * vv * Ip0
                if k > 0 || kp > 0
                    val += 2 * k * kp * vv * Im1
                    val += (kp * d1[b] * v[c] + k * v[b] * d1[c]) * inv1 * (dd * Im1 + Ip1)
                    val += (kp * d2[b] * v[c] + k * v[b] * d2[c]) * inv2 * (-dd * Im1 + Ip1)
                end
                Acc[b, c] += w * val
            end
        end
    end
    A = sum(partial)
    for b = 1:nb, c = 1:(b-1)
        A[b, c] = A[c, b]
    end
    return A
end

# Lift-lift and lift-correction blocks by 3D quadrature, in precision T.
function lift_blocks(lifts, terms, rule3, Rs, E, Z, g, ::Type{T}) where {T}
    nl = length(lifts)
    nb = length(terms)
    nn = length(rule3.w)
    nt = Threads.nthreads()
    Pll = [zeros(T, nl, nl) for _ = 1:nt]
    Plb = [zeros(T, nl, nb) for _ = 1:nt]
    Threads.@threads for chunk = 1:nt
        L0 = zeros(T, nl, 4)
        B0 = zeros(T, nb, 4)
        for node = chunk:nt:nn
            x1, x2, x12, w = rule3.r1[node], rule3.r2[node], rule3.r12[node], rule3.w[node]
            for a = 1:nl
                L0[a, 1], L0[a, 2], L0[a, 3], L0[a, 4] = lifts[a](x1, x2, x12)
            end
            for b = 1:nb
                B0[b, 1], B0[b, 2], B0[b, 3], B0[b, 4] = term_vg_sym(terms[b], x1, x2, x12, T(Rs))
            end
            c1 = (x1^2 + x12^2 - x2^2) / (2x1 * x12)
            c2 = (x2^2 + x12^2 - x1^2) / (2x2 * x12)
            pot = 2 * (-T(Z) / x1 - T(Z) / x2 + T(g) / x12 - T(E))
            form(X, a, Y, b) =
                X[a, 2] * Y[b, 2] + X[a, 3] * Y[b, 3] + 2 * X[a, 4] * Y[b, 4] +
                c1 * (X[a, 2] * Y[b, 4] + X[a, 4] * Y[b, 2]) +
                c2 * (X[a, 3] * Y[b, 4] + X[a, 4] * Y[b, 3]) + pot * X[a, 1] * Y[b, 1]
            for a = 1:nl
                for a2 = a:nl
                    Pll[chunk][a, a2] += w * form(L0, a, L0, a2)
                end
                for b = 1:nb
                    Plb[chunk][a, b] += w * form(L0, a, B0, b)
                end
            end
        end
    end
    All = sum(Pll)
    for a = 1:nl, a2 = 1:(a-1)
        All[a, a2] = All[a2, a]
    end
    return All, sum(Plb)
end

# Correction terms evaluated in T (term_vg in steklov.jl is generic except for the Float64
# fields of Term, which promote correctly).

function exterior_steklov_dd(Att, Atb, Abb; tol = 1e-28)
    T = eltype(Abb)
    sc = 1 ./ sqrt.(abs.(diag(Abb)))
    As = Symmetric(sc .* Abb .* sc')
    ch = cholesky(As; check = false)
    if issuccess(ch)
        X = ch.L \ (sc .* Atb')
        S = Att - X' * X
        info = (; method = :cholesky, kept = size(Abb, 1), total = size(Abb, 1))
    else
        e = eigen(As)
        keep = e.values .> T(tol) * maximum(e.values)
        Y = e.vectors[:, keep]' * (sc .* Atb')
        S = Att - Y' * Diagonal(1 ./ e.values[keep]) * Y
        info = (; method = :eigen, kept = count(keep), total = length(keep))
    end
    return (; S = (S + S') / 2, info...)
end

# Channel-type boundary-layer corrections, vanishing on the sphere:
#   exp(S) Y_i(Omega) x^(p+1) (Rs/R)^(K_i + 4),  x = (R - Rs)/R,  p = 0, ..., P-1.
# High-K trace components decay like (Rs/R)^(K+4) away from the sphere; Hylleraas
# functions with O(1) exponents cannot follow that quickly.
function channel_value(basis, i, p, x, Rs, Z, g)
    r1, r2, r12 = x
    R2 = r1^2 + r2^2
    R = sqrt(R2)
    q = (r1^2 - r2^2) / R2
    s = 1 - r12^2 / R2
    Y = harmonic_value(basis, i, q, s)
    S = -Z * (r1 + r2) + g / 2 * r12
    return exp(S) * Y * ((R - Rs) / R)^(p + 1) * (Rs / R)^(basis.K[i] + 4)
end

function channel_vg(basis, i, p, r1, r2, r12, Rs, Z, g)
    f(x) = channel_value(basis, i, p, x, Rs, Z, g)
    x = [r1, r2, r12]
    grad = ForwardDiff.gradient(f, x)
    return f(x), grad[1], grad[2], grad[3]
end

# Exterior Steklov form on the trace space T_L.
#   zsets:  Hylleraas sets (za, zb, omega), integrated analytically in r12
#   nchan:  P, radial count of channel corrections; lchan: their harmonic degree (default L)
function exterior_form_dd(model, E, Rs, nL; zsets, nchan = 0, lchan = nothing,
                          lift = :gauge, gamma = 1.5, nR = 40, nphi = 24, nu = 20,
                          beta = 2.5, edges = [0.25, 0.75, 2.0, 5.0], npanel = 16, T = Double64,
                          verbose = false, fast = true)
    Z, g = model.Z, model.g
    terms = Term[]
    for (za, zb, om) in zsets
        append!(terms, hylleraas_terms(om, za, zb))
    end
    deg = harmonic_degree(model.basis)
    nc = lchan === nothing ? nL : count(deg .<= lchan)
    rule2 = exterior_rule_2d(Rs, T; nR, nphi, beta, edges, npanel)
    rule3 = exterior_rule_3d(Rs, T; nR, nphi, nu, beta, edges, npanel)
    funcs3 = Any[
        (r1, r2, r12) -> lift_vg(model.basis, i, r1, r2, r12, T(Rs), lift, T(Z), T(g),
                                 model.gauge, T(gamma)) for i = 1:nL
    ]
    for i = 1:nc, p = 0:(nchan-1)
        push!(funcs3, (r1, r2, r12) -> channel_vg(model.basis, i, p, r1, r2, r12, T(Rs), T(Z), T(g)))
    end
    n3 = length(funcs3)
    t1 = @elapsed Ahh = correction_block(terms, rule2, Rs, E, Z, g, T)
    if lift == :gauge && fast
        t2 = @elapsed A33, A3h = blocks3d_fast(model.basis, nL, nc, nchan, terms, rule3, Rs,
                                               E, Z, g, T)
    else
        t2 = @elapsed A33, A3h = lift_blocks(funcs3, terms, rule3, Rs, E, Z, g, T)
    end
    # Correction block: [channel corrections; Hylleraas].
    Att = A33[1:nL, 1:nL]
    Atb = hcat(A33[1:nL, (nL+1):end], A3h[1:nL, :])
    Abb = [A33[(nL+1):end, (nL+1):end] A3h[(nL+1):end, :]; A3h[(nL+1):end, :]' Ahh]
    t3 = @elapsed ext = exterior_steklov_dd(Att, Atb, Abb)
    verbose && @printf("  nb=%d (hyl %d, chan %d): assembly %.1fs + %.1fs, schur %.1fs (%s %d/%d)\n",
                       size(Abb, 1), length(terms), n3 - nL, t1, t2, t3, ext.method, ext.kept, ext.total)
    return (; ext..., nb = size(Abb, 1))
end

function run_match_dd(model, E, Rs, L; zsets, nchan = 0, lchan = nothing,
                      N = max(44, maximum(model.basis.K) + 28), setup = nothing, kwargs...)
    if setup === nothing
        setup = interior_setup(model, E, Rs, L; N)
    end
    ext = exterior_form_dd(model, E, Rs, setup.nL; zsets, nchan, lchan, kwargs...)
    res = match(setup, Float64.(ext.S))
    return (; res..., setup, ext, nterms = ext.nb)
end

# ---------------------------------------------------------------------------
# Faster 3D blocks: analytic gradients, shared per-node tables, chunked matmuls.
# ---------------------------------------------------------------------------

# Y, dY/dq, dY/ds for harmonic i, by a two-partial dual evaluation.
function harmonic_qs(basis, i, q::T, s::T) where {T}
    dq = ForwardDiff.Dual{Nothing}(q, one(T), zero(T))
    ds = ForwardDiff.Dual{Nothing}(s, zero(T), one(T))
    y = harmonic_value(basis, i, dq, ds)
    p = ForwardDiff.partials(y)
    return ForwardDiff.value(y), p[1], p[2]
end

# Fill row `a` of V (value, d1, d2, d12 in V[a, col, 1:4]) for the 3D functions:
# lifts exp(S) Y_i (i <= nL), then channel functions exp(S) Y_i x^(p+1) (Rs/R)^(K_i+4).
function eval3d!(V, a, basis, nL, nc, P, x1::T, x2::T, x12::T, Rs, Z, g) where {T}
    R2 = x1^2 + x2^2
    R = sqrt(R2)
    q = (x1^2 - x2^2) / R2
    s = 1 - x12^2 / R2
    eS = exp(-T(Z) * (x1 + x2) + T(g) / 2 * x12)
    dS = (-T(Z), -T(Z), T(g) / 2)
    dq = (4x1 * x2^2 / R2^2, -4x2 * x1^2 / R2^2, zero(T))
    ds = (2x12^2 * x1 / R2^2, 2x12^2 * x2 / R2^2, -2x12 / R2)
    dR = (x1 / R, x2 / R, zero(T))
    xr = (R - T(Rs)) / R
    dxr = T(Rs) / R2                     # d/dR of (R - Rs)/R
    col = 0
    nmax = max(nL, nc)
    for i = 1:nmax
        Y, Yq, Ys = harmonic_qs(basis, i, q, s)
        f = eS * Y
        g1 = eS * (Y * dS[1] + Yq * dq[1] + Ys * ds[1])
        g2 = eS * (Y * dS[2] + Yq * dq[2] + Ys * ds[2])
        g3 = eS * (Y * dS[3] + Yq * dq[3] + Ys * ds[3])
        if i <= nL
            V[a, i, 1], V[a, i, 2], V[a, i, 3], V[a, i, 4] = f, g1, g2, g3
        end
        if i <= nc
            K = basis.K[i]
            base = (T(Rs) / R)^(K + 4)
            dbase = -(K + 4) / R * base
            for p = 0:(P-1)
                rho = xr^(p + 1) * base
                drho = (p + 1) * xr^p * dxr * base + xr^(p + 1) * dbase
                c = nL + (i - 1) * P + p + 1
                V[a, c, 1] = f * rho
                V[a, c, 2] = g1 * rho + f * drho * dR[1]
                V[a, c, 3] = g2 * rho + f * drho * dR[2]
                V[a, c, 4] = g3 * rho
            end
        end
    end
    return nothing
end

function hylleraas_eval!(H, a, terms, x1::T, x2::T, x12::T, Rs) where {T}
    for (b, t) in enumerate(terms)
        H[a, b, 1], H[a, b, 2], H[a, b, 3], H[a, b, 4] = term_vg_sym(t, x1, x2, x12, T(Rs))
    end
    return nothing
end

function blocks3d_fast(basis, nL, nc, P, terms, rule3, Rs, E, Z, g, ::Type{T};
                       chunk = 128) where {T}
    n3 = nL + nc * P
    nb = length(terms)
    nn = length(rule3.w)
    starts = collect(1:chunk:nn)
    nt = Threads.nthreads()
    P33 = [zeros(T, n3, n3) for _ = 1:nt]
    P3h = [zeros(T, n3, nb) for _ = 1:nt]
    Threads.@threads for tid = 1:nt
        V = zeros(T, chunk, n3, 4)
        H = zeros(T, chunk, nb, 4)
        for st in starts[tid:nt:end]
            idx = st:min(st + chunk - 1, nn)
            nc_ = length(idx)
            fill!(V, zero(T))
            fill!(H, zero(T))
            wv = zeros(T, chunk)
            c1 = zeros(T, chunk)
            c2 = zeros(T, chunk)
            pot = zeros(T, chunk)
            for (a, node) in enumerate(idx)
                x1, x2, x12 = rule3.r1[node], rule3.r2[node], rule3.r12[node]
                eval3d!(V, a, basis, nL, nc, P, x1, x2, x12, Rs, Z, g)
                hylleraas_eval!(H, a, terms, x1, x2, x12, Rs)
                wv[a] = rule3.w[node]
                c1[a] = (x1^2 + x12^2 - x2^2) / (2x1 * x12)
                c2[a] = (x2^2 + x12^2 - x1^2) / (2x2 * x12)
                pot[a] = 2 * (-T(Z) / x1 - T(Z) / x2 + T(g) / x12 - T(E))
            end
            V0, V1, V2, V12 = (view(V, :, :, c) for c = 1:4)
            H0, H1, H2, H12 = (view(H, :, :, c) for c = 1:4)
            # form(X, Y) = X1'Y1 + X2'Y2 + 2 X12'Y12 + c1 (X1'Y12 + X12'Y1)
            #              + c2 (X2'Y12 + X12'Y2) + pot X0'Y0, all with weights.
            W1 = wv .* V1 .+ (wv .* c1) .* V12
            W2 = wv .* V2 .+ (wv .* c2) .* V12
            W12 = 2 .* wv .* V12 .+ (wv .* c1) .* V1 .+ (wv .* c2) .* V2
            W0 = (wv .* pot) .* V0
            P33[tid] .+= W1' * V1 .+ W2' * V2 .+ W12' * V12 .+ W0' * V0
            P3h[tid] .+= W1' * H1 .+ W2' * H2 .+ W12' * H12 .+ W0' * H0
        end
    end
    A33 = sum(P33)
    return (A33 + A33') / 2, sum(P3h)
end
