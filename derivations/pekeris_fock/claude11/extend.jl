# Real-analytic continuation of the Fock hierarchy psi_kp across the collinear edge X = -pi/2 (side = -1) or X = +pi/2
# (side = +1) of the (X, Y) = (alpha-beta, alpha+beta) square, up to the next degenerate line X0 = -3pi/2 (resp. +3pi/2),
# where the continued branch may carry logarithms. The operator is separable, w*(Lambda^2 - c) ~ M_X - M_Y with
# M = d(cos d) + (c/16) cos, so in the L_Y eigenbasis each Y-mode is an ODE in X:
#   element 2: [Xmid, X1], plain Chebyshev-Lobatto; analytic at X1, value at X1 from the physical spectral solution;
#   element 1: [X0, Xmid], psi = sum_m A_m(X) l^m, l = log(-cos X), A_m analytic at X0.
# Output: the log coefficients A_m(X0, Y), m >= 1, on the physical Y grid, for each requested psi_kp and Z-sector.
include(joinpath(@__DIR__, "..", "claude9", "fock_spectral_copy.jl"))
using Serialization

lobatto(T, N) = [cos(T(π) * j / (N - 1)) for j = 0:N-1]          # descending from 1 to -1
function lobatto_D(x::Vector{T}) where {T}
    N = length(x)
    c = [(j == 1 || j == N ? T(2) : T(1)) * (-1)^(j - 1) for j = 1:N]
    D = zeros(T, N, N)
    for i = 1:N, j = 1:N
        i == j || (D[i, j] = c[i] / c[j] / (x[i] - x[j]))
    end
    for i = 1:N
        D[i, i] = -sum(D[i, :])
    end
    D
end

struct Elem{T}
    X::Vector{T}      # physical X nodes
    D::Matrix{T}      # d/dX
    iA::Int           # index of the node at the "A" end (X0 for element 1, X1 for element 2)
    iB::Int           # index of the node at Xmid
end
function Elem(T, a, b, N)   # nodes from a (x = -1 ... ) to b
    x = lobatto(T, N)                    # x[1] = 1, x[N] = -1
    X = (a + b) / 2 .+ (b - a) / 2 .* x  # X[1] = b, X[N] = a
    D = lobatto_D(x) .* (2 / (b - a))
    Elem{T}(X, D, N, 1)                  # a end is node N, b end is node 1
end

function run_extend(; T = BigFloat, n = 80, N = 60, Zs = 0:4, kmax = 4, side = -1, E = 0)
    P = T(π)
    g = SGrid(T, n)
    X1 = side * P / 2; X0 = side * 3P / 2; Xm = side * P
    # element 1 runs between X0 and Xm, element 2 between Xm and X1; A end = X0 / X1, B end = Xm
    e1 = Elem(T, X0, Xm, N); e2 = Elem(T, X1, Xm, N)
    Ys = g.Y
    # physical interpolation weights to X = X1 (tx = side)
    xn = [cos(P * (j + T(1) / 2) / n) for j = 0:n-1]
    bw = [(-1)^j * sin(P * (j + T(1) / 2) / n) for j = 0:n-1]
    ax = bw ./ (T(side) .- xn); ax ./= sum(ax)
    grid(e) = ([(x + y) / 2 for x in e.X, y in Ys], [(y - x) / 2 for x in e.X, y in Ys])
    results = Dict{Any,Any}()
    for Z in Zs
        Zt = T(Z)
        wVf(α, β) = 2sqrt(T(2)) * sin(α) * cos(β / 2) - 4Zt * sin(β) * (cos(α / 2) + sin(α / 2))
        wf(α, β) = 2sin(α) * sin(β)
        wVp = wV_phys(g, Zt)
        a1, b1 = grid(e1); a2, b2 = grid(e2)
        w1 = wf.(a1, b1); wV1 = wVf.(a1, b1); w2 = wf.(a2, b2); wV2 = wVf.(a2, b2)
        ψ = Dict{Tuple{Int,Int},Matrix{T}}((0, 0) => ones(T, n, n))
        # extended: ψ2[(k,p)] :: N x n values on element 2; ψ1[(k,p)] :: Vector of N x n (coefficient of l^m, m = 0..)
        ψ2 = Dict{Tuple{Int,Int},Matrix{T}}((0, 0) => ones(T, N, n))
        ψ1 = Dict{Tuple{Int,Int},Vector{Matrix{T}}}((0, 0) => [ones(T, N, n)])
        get0(k, p) = get(ψ, (k, p), zeros(T, n, n))
        get2(k, p) = get(ψ2, (k, p), zeros(T, N, n))
        get1(k, p) = get(ψ1, (k, p), [zeros(T, N, n)])
        padd(a, b) = [(m <= length(a) ? a[m] : zero(a[1])) .+ (m <= length(b) ? b[m] : zero(b[1])) for m = 1:max(length(a), length(b))]
        pscale(s, a) = [s .* x for x in a]
        for k = 1:kmax
            c = k * (k + 4)
            o = SOps(g, c)
            H = iseven(k) ? sharmonics(g, k) : Matrix{T}[]
            # harmonics on the element grids
            function harm_on(α, β)
                x1 = cos.(α); x2 = cos.(β); s2 = 1 .- x1 .^ 2
                Hs = Matrix{T}[]
                for l = 0:k÷2
                    m = k ÷ 2 - l
                    C0 = ones(T, size(x1)); C1 = 2(l + 1) * x1; Cm = m == 0 ? C0 : C1
                    for j = 2:m
                        Cn = (2(j + l) * x1 .* C1 - (j + 2l) * C0) / j; C0, C1 = C1, Cn; Cm = Cn
                    end
                    Pl = zeros(T, size(x1))
                    for kk = 0:l÷2
                        Pl .+= T((-1)^kk * binomial(l, kk) * binomial(2l - 2kk, l)) / 2^l .* x2 .^ (l - 2kk) .* s2 .^ kk
                    end
                    push!(Hs, Cm .* Pl)
                end
                Hs
            end
            H1 = iseven(k) ? harm_on(a1, b1) : Matrix{T}[]
            H2 = iseven(k) ? harm_on(a2, b2) : Matrix{T}[]
            # per-mode LU factorizations for M_X - λ on each element
            MX(e) = Diagonal(cos.(e.X)) * e.D^2 - Diagonal(sin.(e.X)) * e.D + (T(c) / 16) * Diagonal(cos.(e.X))
            M1 = MX(e1); M2 = MX(e2)
            function sysmat(M, e, λ)
                A = M - λ * I
                A = Matrix(A)
                A[e.iB, :] .= 0; A[e.iB, e.iA] = 1      # replace the Xmid row by "value at the A end"
                lu(A)
            end
            LU1 = [sysmat(M1, e1, o.λY[j]) for j = 1:n]
            LU2 = [sysmat(M2, e2, o.λY[j]) for j = 1:n]
            sX1 = sin.(e1.X); cX1 = cos.(e1.X)
            ℓm = log(-cos(Xm)); tm = tan(Xm)
            for p = (k÷2):-1:0
                wh(ψk1) = g.w .* (2(k + 2) * (p + 1) * ψk1 + (p + 1) * (p + 2) * get0(k, p + 2) + 2E * get0(k - 2, p)) - 2wVp .* get0(k - 1, p)
                wh2(f) = w2 .* (2(k + 2) * (p + 1) * f + (p + 1) * (p + 2) * get2(k, p + 2) + 2E * get2(k - 2, p)) - 2wV2 .* get2(k - 1, p)
                wh1(f) = padd([w1 .* x for x in padd(padd(pscale(2(k + 2) * (p + 1), f), pscale((p + 1) * (p + 2), get1(k, p + 2))), pscale(2E, get1(k - 2, p)))],
                    pscale(-2, [wV1 .* y for y in get1(k - 1, p)]))
                if iseven(k) && p < k ÷ 2 && haskey(ψ, (k, p + 1))
                    G = [sip(g, a, b) for a in H, b in H]
                    s = [sum(g.W ./ g.w .* h .* wh(get0(k, p + 1))) for h in H]
                    cc = -(G \ s) ./ (2(k + 2) * (p + 1))
                    ψ[(k, p + 1)] += sum(cc[l] * H[l] for l in eachindex(H))
                    ψ2[(k, p + 1)] += sum(cc[l] * H2[l] for l in eachindex(H))
                    v = ψ1[(k, p + 1)]; v[1] = v[1] + sum(cc[l] * H1[l] for l in eachindex(H))
                end
                ψkp = ssolve(g, o, wh(get0(k, p + 1)), H)
                ψ[(k, p)] = ψkp
                # ---- extension ----
                r2 = (-wh2(get2(k, p + 1)) ./ 16) * o.iVY'                 # N x n modal
                r1 = [(-x ./ 16) * o.iVY' for x in wh1(get1(k, p + 1))]   # list over powers m
                v1 = vec((ax' * ψkp) * o.iVY')                             # modal values at X1
                K = length(r1)   # source has powers 0..K-1; solution may have 0..K
                Am = [zeros(T, N, n) for _ = 0:K]
                S2 = zeros(T, N, n)
                for j = 1:n
                    # element 2
                    rhs = copy(r2[:, j]); rhs[e2.iB] = v1[j]
                    y2 = LU2[j] \ rhs
                    S2[:, j] = y2
                    f_m = y2[e2.iB]; d_m = (e2.D*y2)[e2.iB]
                    # element 1, top-down
                    A = [zeros(T, N) for _ = 0:K]
                    function rhs1(mm, extra)
                        r = mm + 1 <= K ? copy(r1[mm+1][:, j]) : zeros(T, N)
                        r .-= extra
                        if mm + 1 <= K
                            Ap = A[mm+2]; dAp = e1.D * Ap
                            r .-= (mm + 1) .* (2 .* sX1 .* dAp .+ cX1 .* Ap)
                        end
                        if mm + 2 <= K
                            App = A[mm+3]; dApp = e1.D * App
                            q = similar(App)
                            for i = 1:N
                                q[i] = i == e1.iA ? -sX1[i] * dApp[i] : sX1[i]^2 * App[i] / cX1[i]
                            end
                            r .-= (mm + 2) * (mm + 1) .* q
                        end
                        r
                    end
                    for mm = K:-1:2
                        r = rhs1(mm, zeros(T, N)); r[e1.iB] = 0
                        A[mm+1] = LU1[j] \ r
                    end
                    r = rhs1(1, zeros(T, N)); r[e1.iB] = 0; P1 = LU1[j] \ r
                    h = zeros(T, N); h[e1.iB] = 1; φ = LU1[j] \ h
                    A[2] = P1
                    r = rhs1(0, zeros(T, N)); r[e1.iB] = 0; P0 = LU1[j] \ r
                    dφ = e1.D * φ
                    r = -(2 .* sX1 .* dφ .+ cX1 .* φ); r[e1.iB] = 0; χ = LU1[j] \ r
                    # psi at Xm: A0 + A1 l + sum_{m>=2} Am l^m ; derivative
                    val(Avec) = sum(Avec[m+1][e1.iB] * ℓm^m for m = 0:K)
                    der(Avec) = sum((e1.D*Avec[m+1])[e1.iB] * ℓm^m + (m > 0 ? m * Avec[m+1][e1.iB] * ℓm^(m - 1) * tm : zero(T)) for m = 0:K)
                    Abase = deepcopy(A); Abase[1] = P0
                    Ac1 = [zeros(T, N) for _ = 0:K]; Ac1[2] = φ; Ac1[1] = χ
                    Ac0 = [zeros(T, N) for _ = 0:K]; Ac0[1] = φ
                    Mm = [val(Ac1) val(Ac0); der(Ac1) der(Ac0)]
                    cvec = Mm \ [f_m - val(Abase), d_m - der(Abase)]
                    for m = 0:K
                        Am[m+1][:, j] = Abase[m+1] + cvec[1] * Ac1[m+1] + cvec[2] * Ac0[m+1]
                    end
                end
                ψ2[(k, p)] = S2 * o.VY'
                # drop trailing zero powers
                vals = [x * o.VY' for x in Am]
                while length(vals) > 1 && maximum(abs.(vals[end])) < T(10)^(-precision(T) * 0.25)
                    pop!(vals)
                end
                ψ1[(k, p)] = vals
                println("Z=$Z (k,p)=($k,$p) log degree on element 1: ", length(vals) - 1,
                    "  max|A_m(X0)|: ", [Float64(maximum(abs.(v[e1.iA, :]))) for v in vals])
                flush(stdout)
            end
        end
        results[Z] = (ψ = ψ, ψ1 = ψ1, ψ2 = ψ2)
    end
    (g = g, e1 = e1, e2 = e2, results = results, side = side)
end
wV_phys(g::SGrid{T}, Z) where {T} = 2sqrt(T(2)) * sin.(g.α) .* cos.(g.β ./ 2) - 4Z * sin.(g.β) .* (cos.(g.α ./ 2) + sin.(g.α ./ 2))

if abspath(PROGRAM_FILE) == @__FILE__
    bits = parse(Int, get(ENV, "BITS", "256")); setprecision(BigFloat, bits)
    n = parse(Int, get(ENV, "NGRID", "80")); N = parse(Int, get(ENV, "NEL", "60"))
    side = parse(Int, get(ENV, "SIDE", "-1")); kmax = parse(Int, get(ENV, "KMAX", "4"))
    Zs = eval(Meta.parse(get(ENV, "ZS", "0:4")))
    t = @elapsed R = run_extend(; n, N, side, kmax, Zs)
    println("elapsed ", t)
    serialize(get(ENV, "OUT", joinpath(@__DIR__, "ext_side$(side)_n$(n)_N$(N).jls")), R)
end
