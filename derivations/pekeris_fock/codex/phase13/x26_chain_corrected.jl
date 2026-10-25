# Iterated circle variations of the continued Fock hierarchy on the strip J x Y_phys, J = [Xa, Xb] = side*[pi/2, 3pi/2]
# (between the physical edge Xa and the next degenerate line Xb). V_0 = continued psi (analytic at Xa); V_{s+1} = N V_s,
# the log coefficient (d/dl of the l-polynomial, l = log(-cos X)) at the far end of step s, which is analytic there and
# is continued back across J. N commutes with the hierarchy operator and with multiplication by the analytic w, wV, so
# every V_s obeys the same recursion; only V_0 carries the Y4 harmonic additions. Per Y-mode: element S (start end,
# plain Chebyshev-Lobatto, analytic at the start end with prescribed value), element F (far end, sum_m A_m l^m).
include(joinpath(@__DIR__, "../../claude11/extend.jl"))   # lobatto, Elem, SGrid, SOps, sharmonics, ssolve, sip, wV_phys
using Serialization

function harm_grid(α, β, k, T)
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
padd(a, b) = [(m <= length(a) ? a[m] : zero(a[1])) .+ (m <= length(b) ? b[m] : zero(b[1])) for m = 1:max(length(a), length(b))]
pscale(s, a) = [s .* x for x in a]

function run_chain(; T = BigFloat, n = 80, N = 60, Zs = 0:4, kmax = 4, side = 1, steps = 4, E = 0, wvk = nothing)
    P = T(π)
    g = SGrid(T, n)
    Xa = side * P / 2; Xb = side * 3P / 2; Xm = side * P
    eA = Elem(T, Xa, Xm, N); eB = Elem(T, Xb, Xm, N)   # A end of eA is Xa, of eB is Xb; B end of both is Xm
    Ys = g.Y
    xn = [cos(P * (j + T(1) / 2) / n) for j = 0:n-1]
    bw = [(-1)^j * sin(P * (j + T(1) / 2) / n) for j = 0:n-1]
    axa = bw ./ (T(side) .- xn); axa ./= sum(axa)
    grid(e) = ([(x + y) / 2 for x in e.X, y in Ys], [(y - x) / 2 for x in e.X, y in Ys])
    ops = Dict(k => SOps(g, k * (k + 4)) for k = 1:kmax)
    MX(e, c) = Diagonal(cos.(e.X)) * e.D^2 - Diagonal(sin.(e.X)) * e.D + (T(c) / 16) * Diagonal(cos.(e.X))
    function sysmat(M, e, λ)
        A = Matrix(M - λ * I)
        A[e.iB, :] .= 0; A[e.iB, e.iA] = 1
        lu(A)
    end
    LU = Dict()
    for k = 1:kmax, (nm, e) in ((:A, eA), (:B, eB))
        M = MX(e, k * (k + 4))
        LU[(k, nm)] = [sysmat(M, e, ops[k].λY[j]) for j = 1:n]
    end
    out = Dict{Any,Any}()
    for Z in Zs
        Zt = T(Z)
        wv0p = 2sqrt(T(2)) .* sin.(g.α) .* cos.(g.β ./ 2); wv1p = 4 .* sin.(g.β) .* (cos.(g.α ./ 2) .+ sin.(g.α ./ 2))
        wVk(k) = wvk === nothing ? wv0p .- Zt .* wv1p : wvk(k)[1] .* wv0p .- wvk(k)[2] .* wv1p
        geo = Dict()
        for (nm, e) in ((:A, eA), (:B, eB))
            α, β = grid(e)
            geo[nm] = (w = 2 .* sin.(α) .* sin.(β), wv0 = 2sqrt(T(2)) .* sin.(α) .* cos.(β ./ 2), wv1 = 4 .* sin.(β) .* (cos.(α ./ 2) .+ sin.(α ./ 2)),
                H = Dict(k => harm_grid(α, β, k, T) for k = 2:2:kmax), sX = sin.(e.X), cX = cos.(e.X))
        end
        wVe(gg, k) = wvk === nothing ? gg.wv0 .- Zt .* gg.wv1 : wvk(k)[1] .* gg.wv0 .- wvk(k)[2] .* gg.wv1
        # physical hierarchy (needed for step-0 start values and harmonic constants)
        ψ = Dict{Tuple{Int,Int},Matrix{T}}((0, 0) => ones(T, n, n))
        get0(k, p) = get(ψ, (k, p), zeros(T, n, n))
        startvals = Dict{Tuple{Int,Int},Vector{T}}()   # per order: values at the start end (Y grid)
        harmc = Dict{Tuple{Int,Int},Vector{T}}()       # cc added to (k,p+1)
        for k = 1:kmax
            o = ops[k]
            H = iseven(k) ? sharmonics(g, k) : Matrix{T}[]
            for p = (k÷2):-1:0
                wh(ψk1) = g.w .* (2(k + 2) * (p + 1) * ψk1 + (p + 1) * (p + 2) * get0(k, p + 2) + 2E * get0(k - 2, p)) - 2wVk(k) .* get0(k - 1, p)
                if iseven(k) && p < k ÷ 2 && haskey(ψ, (k, p + 1))
                    G = [sip(g, a, b) for a in H, b in H]
                    s = [sum(g.W ./ g.w .* h .* wh(get0(k, p + 1))) for h in H]
                    cc = -(G \ s) ./ (2(k + 2) * (p + 1))
                    ψ[(k, p + 1)] += sum(cc[l] * H[l] for l in eachindex(H))
                    harmc[(k, p + 1)] = cc
                end
                ψ[(k, p)] = ssolve(g, o, wh(get0(k, p + 1)), H)
            end
        end
        for s = 0:steps
            S, F = iseven(s) ? (:A, :B) : (:B, :A)
            eS = S == :A ? eA : eB; eF = F == :A ? eA : eB
            gS = geo[S]; gF = geo[F]
            VS = Dict{Tuple{Int,Int},Matrix{T}}()                  # element S values
            VF = Dict{Tuple{Int,Int},Vector{Matrix{T}}}()          # element F l-coefficients
            if s == 0
                VS[(0, 0)] = ones(T, N, n); VF[(0, 0)] = [ones(T, N, n)]
            end
            gS0(k, p) = get(VS, (k, p), zeros(T, N, n))
            gF0(k, p) = get(VF, (k, p), [zeros(T, N, n)])
            newstart = Dict{Tuple{Int,Int},Vector{T}}()
            maxdeg = 0
            for k = 1:kmax
                o = ops[k]; c = k * (k + 4)
                ℓm = log(-cos(Xm)); tm = tan(Xm)
                for p = (k÷2):-1:0
                    if false # x26: final physical startvals already contain this correction
                        cc = harmc[(k, p + 1)]
                        VS[(k, p + 1)] += sum(cc[l] * gS.H[k][l] for l in eachindex(cc))
                        v = VF[(k, p + 1)]; v[1] = v[1] + sum(cc[l] * gF.H[k][l] for l in eachindex(cc))
                    end
                    whS(f) = gS.w .* (2(k + 2) * (p + 1) * f + (p + 1) * (p + 2) * gS0(k, p + 2) + 2E * gS0(k - 2, p)) - 2wVe(gS, k) .* gS0(k - 1, p)
                    whF(f) = padd([gF.w .* x for x in padd(padd(pscale(2(k + 2) * (p + 1), f), pscale((p + 1) * (p + 2), gF0(k, p + 2))), pscale(2E, gF0(k - 2, p)))],
                        pscale(-2, [wVe(gF, k) .* y for y in gF0(k - 1, p)]))
                    rS = (-whS(gS0(k, p + 1)) ./ 16) * o.iVY'
                    rF = [(-x ./ 16) * o.iVY' for x in whF(gF0(k, p + 1))]
                    sv = s == 0 ? vec(axa' * ψ[(k, p)]) : get(startvals, (k, p), zeros(T, n))
                    v1 = vec(sv' * o.iVY')
                    K = length(rF)
                    Am = [zeros(T, N, n) for _ = 0:K]
                    SS = zeros(T, N, n)
                    sXF = gF.sX; cXF = gF.cX
                    for j = 1:n
                        rhs = copy(rS[:, j]); rhs[eS.iB] = v1[j]
                        y2 = LU[(k, S)][j] \ rhs
                        SS[:, j] = y2
                        f_m = y2[eS.iB]; d_m = (eS.D*y2)[eS.iB]
                        A = [zeros(T, N) for _ = 0:K]
                        function rhs1(mm)
                            r = mm + 1 <= K ? copy(rF[mm+1][:, j]) : zeros(T, N)
                            if mm + 1 <= K
                                Ap = A[mm+2]; dAp = eF.D * Ap
                                r .-= (mm + 1) .* (2 .* sXF .* dAp .+ cXF .* Ap)
                            end
                            if mm + 2 <= K
                                App = A[mm+3]; dApp = eF.D * App
                                q = [i == eF.iA ? -sXF[i] * dApp[i] : sXF[i]^2 * App[i] / cXF[i] for i = 1:N]
                                r .-= (mm + 2) * (mm + 1) .* q
                            end
                            r
                        end
                        LUF = LU[(k, F)][j]
                        for mm = K:-1:2
                            r = rhs1(mm); r[eF.iB] = 0
                            A[mm+1] = LUF \ r
                        end
                        r = rhs1(1); r[eF.iB] = 0; A[2] = LUF \ r
                        h = zeros(T, N); h[eF.iB] = 1; φ = LUF \ h
                        r = rhs1(0); r[eF.iB] = 0; P0 = LUF \ r
                        dφ = eF.D * φ
                        r = -(2 .* sXF .* dφ .+ cXF .* φ); r[eF.iB] = 0; χ = LUF \ r
                        val(Av) = sum(Av[m+1][eF.iB] * ℓm^m for m = 0:K)
                        der(Av) = sum((eF.D*Av[m+1])[eF.iB] * ℓm^m + (m > 0 ? m * Av[m+1][eF.iB] * ℓm^(m - 1) * tm : zero(T)) for m = 0:K)
                        Ab = deepcopy(A); Ab[1] = P0
                        Ac1 = [zeros(T, N) for _ = 0:K]; Ac1[2] = φ; Ac1[1] = χ
                        Ac0 = [zeros(T, N) for _ = 0:K]; Ac0[1] = φ
                        cv = [val(Ac1) val(Ac0); der(Ac1) der(Ac0)] \ [f_m - val(Ab), d_m - der(Ab)]
                        for m = 0:K
                            Am[m+1][:, j] = Ab[m+1] + cv[1] * Ac1[m+1] + cv[2] * Ac0[m+1]
                        end
                    end
                    VS[(k, p)] = SS * o.VY'
                    vals = [x * o.VY' for x in Am]
                    while length(vals) > 1 && maximum(abs.(vals[end])) < T(10)^(-50)
                        pop!(vals)
                    end
                    VF[(k, p)] = vals
                    deg = length(vals) - 1
                    maxdeg = max(maxdeg, deg)
                    # N at the far end: d/dl; start value of the next step = coefficient of l^1 at the far node
                    newstart[(k, p)] = deg >= 1 ? vals[2][eF.iA, :] : zeros(T, n)
                    if deg >= 2
                        println("  WARNING Z=$Z s=$s ($k,$p): log degree $deg at far end; N V not analytic there")
                    end
                end
            end
            out[(Z, s)] = newstart
            println("Z=$Z step $s (start $(S), far $(F)): max log degree $maxdeg; |N V| per order: ",
                join(["$(kp)=>" * string(round(Float64(maximum(abs.(v))), sigdigits = 4)) for (kp, v) in sort(collect(newstart))], " "))
            flush(stdout)
            startvals = newstart
        end
    end
    (g = g, out = out, side = side, steps = steps)
end

if abspath(PROGRAM_FILE) == @__FILE__
    setprecision(BigFloat, parse(Int, get(ENV, "BITS", "256")))
    n = parse(Int, get(ENV, "NGRID", "80")); N = parse(Int, get(ENV, "NEL", "60"))
    steps = parse(Int, get(ENV, "STEPS", "6")); Zs = eval(Meta.parse(get(ENV, "ZS", "0:4")))
    lad = get(ENV, "LADDER", "")
    wvk = lad == "" ? nothing : (k -> eval(Meta.parse(lad))[k])
    t = @elapsed R = run_chain(; n, N, steps, Zs, side = parse(Int, get(ENV, "SIDE", "1")), wvk)
    println("elapsed ", t)
    serialize(get(ENV, "OUT", joinpath(@__DIR__, "chain_n$(n)_N$(N).jls")), R)
end
if abspath(PROGRAM_FILE) == @__FILE__
    println("x26 peak_rss_bytes ", Sys.maxrss())
end
