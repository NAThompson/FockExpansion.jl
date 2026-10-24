# Exterior log-derivative of the cusp-gauge channel ODE as an inverse factorial series.
#
#   Phi'' + (5/R + 2U) Phi' + (-Lambda/R^2 + D/R + T) Phi = 0,  T = W + 2E.
#
# Y = Phi' Phi^{-1} for the span of the m subdominant solutions solves the Riccati equation
#   Y' + Y^2 + (5/R + 2U) Y - Lambda/R^2 + D/R + T = 0
# with formal series Y = sum_n Y_n R^{-n}:
#   Y0^2 + 2U Y0 + T = 0                                   (slow solvent)
#   (Y0 + 2U) Y1 + Y1 Y0 = -(5 Y0 + D)
#   (Y0 + 2U) Yn + Yn Y0 = (n - 6) Y_{n-1} - sum_{j=1}^{n-1} Y_j Y_{n-j} + Lambda [n == 2].
# Borel: Y - Y0 = int_0^inf exp(-R t) yhat(t) dt, yhat = sum_{n>=1} Y_n t^{n-1}/(n-1)!.
# With t = -log(1-u)/omega and L^j/j! = sum_k |s(k,j)| u^k/k!,
#   Y(R) = Y0 + sum_k B_k / (x)_{k+1},  x = R/omega,
#   B_k = sum_{j=0}^k |s(k,j)| omega^{-j-1} Y_{j+1}.
include("/Users/nathompson7/Documents/Codex/2026-09-05/how-much-compu/outputs/boundary_matching/cusp_matching.jl")

sylvester_operator(A, B) = kron(Matrix{eltype(A)}(I, size(B)...), A) + kron(transpose(B), Matrix{eltype(A)}(I, size(A)...))

"""Slow solvent Y0 in BigFloat by Newton refinement of cusptail's Float64 solvent."""
function slow_solvent(b, E; prec = 512)
    setprecision(BigFloat, prec) do
        Y, _, _ = cusptail(b, E)
        U = big.(b.U)
        T = big.(b.W) + 2big(E) * I
        Y = big.(Y)
        m = size(Y, 1)
        for _ = 1:12
            Q = Y * Y + 2U * Y + T
            norm(Q) < big(2.0)^(-prec + 20) && break
            S = sylvester_operator(Y + 2U, Y)
            Y -= reshape(S \ vec(Q), m, m)
        end
        Y
    end
end

"""Formal Riccati coefficients Y_0..Y_N (BigFloat)."""
function riccati_coeffs(b, E, N; prec = 512)
    setprecision(BigFloat, prec) do
        Y0 = slow_solvent(b, E; prec)
        U = big.(b.U)
        D = big.(b.D)
        Lam = Diagonal(big.(b.lam))
        m = size(Y0, 1)
        S = lu(sylvester_operator(Y0 + 2U, Y0))
        Ys = [Y0]
        push!(Ys, reshape(S \ vec(-(5Y0 + D)), m, m))
        for n = 2:N
            C = (n - 6) * Ys[n]
            for j = 1:n-1
                C -= Ys[j+1] * Ys[n-j+1]
            end
            n == 2 && (C += Lam)
            push!(Ys, reshape(S \ vec(C), m, m))
        end
        Ys
    end
end

"""Unsigned Stirling numbers of the first kind c[k+1, j+1] = |s(k, j)|, 0 <= j <= k <= N."""
function stirling1(N)
    c = zeros(BigInt, N + 1, N + 1)
    c[1, 1] = 1
    for k = 0:N-1, j = 0:k+1
        c[k+2, j+1] = k * c[k+1, j+1] + (j >= 1 ? c[k+1, j] : 0)
    end
    c
end

"""Inverse factorial series coefficients B_0..B_{K} (needs Ys[1..K+2])."""
function factorial_coeffs(Ys, omega, K; prec = precision(Ys[1][1]), stir = stirling1(K))
    setprecision(BigFloat, prec) do
        w = big(omega)
        Bk = Vector{Matrix{BigFloat}}(undef, K + 1)
        for k = 0:K
            acc = zero(Ys[1])
            for j = 0:k
                acc += (stir[k+1, j+1] * w^(-j - 1)) * Ys[j+2]
            end
            Bk[k+1] = acc
        end
        Bk
    end
end

"""Partial sums Y_K(R) for K = 0..length(Bk)-1, returned as Float64 matrices."""
function factorial_partial_sums(Y0, Bk, omega, R; prec = precision(Y0[1]))
    setprecision(BigFloat, prec) do
        x = big(R) / big(omega)
        poch = x
        S = copy(Y0)
        sums = Matrix{Float64}[]
        for k = 0:length(Bk)-1
            S += Bk[k+1] / poch
            push!(sums, Float64.(S))
            poch *= x + k + 1
        end
        sums
    end
end

"""Truncated asymptotic sum sum_{n<=K} Y_n R^{-n} partial sums."""
function asymptotic_partial_sums(Ys, R)
    S = zero(Ys[1])
    out = Matrix{Float64}[]
    for n = 0:length(Ys)-1
        S += Ys[n+1] / big(R)^n
        push!(out, Float64.(S))
    end
    out
end

"""Gap between slow and fast pencil roots (Borel radius of the Y series)."""
function pencil_roots(b, E)
    m = length(b.K)
    T = b.W + 2E * I
    ev = eigvals([zeros(m, m) I; -T -2b.U])
    ev = ev[sortperm(real(ev))]
    ev[1:m], ev[m+1:end]
end

"""Formal series of the Moebius chart W = (cI - Y)^{-1} from the Y series."""
function chart_coeffs(Ys, c)
    m = size(Ys[1], 1)
    A0inv = inv(c * I - Ys[1])
    Ws = [A0inv]
    for n = 1:length(Ys)-1
        acc = zero(Ys[1])
        for j = 1:n
            acc += Ys[j+1] * Ws[n-j+1]
        end
        push!(Ws, A0inv * acc)
    end
    Ws
end

"""Inverse factorial series value of Y at R using K+1 terms (Float64 result)."""
function factorial_Y(Ys, Bk, omega, R)
    factorial_partial_sums(Ys[1], Bk, omega, R)[end]
end

"""
One Taylor expansion of the fundamental matrix about R1 (radius of convergence R1),
started from the boundary map Y(R1) and evaluated at R0. Same recurrence as one
cuspexterior panel, but a single expansion of high order in BigFloat.
Returns the partial-sum boundary maps for orders 1..order.
"""
function taylor_connect(b, E, Y1, R1, R0, order; prec = 256)
    setprecision(BigFloat, prec) do
        m = length(b.K)
        r = big(R1)
        eye = Matrix{BigFloat}(I, m, m)
        U = big.(b.U); Dm = big.(b.D); T = big.(b.W) + 2big(E) * I; lam = Diagonal(big.(b.lam))
        B = [zeros(BigFloat, m, m) for _ = 0:order]
        B[1] = eye; B[2] = r * big.(Y1)
        for n = 0:order-2
            rhs = ((n + 1) * (2n + 5) * eye + 2r * (n + 1) * U) * B[n+2]
            rhs += (n * (n + 4) * eye + 4r * n * U - lam + r * Dm + r^2 * T) * B[n+1]
            n >= 1 && (rhs += (2r * (n - 1) * U + r * Dm + 2r^2 * T) * B[n])
            n >= 2 && (rhs += r^2 * T * B[n-1])
            B[n+3] = -rhs / ((n + 1) * (n + 2))
        end
        x = (big(R0) - r) / r
        P = copy(eye); Q = zeros(BigFloat, m, m); out = Matrix{Float64}[]
        for n = 1:order
            P += B[n+1] * x^n
            Q += n * B[n+1] * x^(n - 1)
            push!(out, Float64.((Q / r) / P))
        end
        out
    end
end

"""Scan for singularities of Y in Re R > rmin, |Im R| <= ymax by complex Riccati integration."""
function pole_scan(b, E, Ys; ymax = 10.0, dy = 0.02, rstart = 14.0, rmin = 0.3, h = 0.004)
    yasym(R) = sum(ComplexF64.(Ys[n+1]) / R^n for n = 0:min(40, length(Ys) - 1))
    U = b.U; Dm = b.D; T = b.W + 2E * I; Lam = Diagonal(b.lam)
    f(R, Y) = -(Y * Y + (5 / R * I + 2U) * Y - Lam / R^2 + Dm / R + T)
    # per line: (Im R, max |Y|, Re R at the max)
    lines = Tuple{Float64,Float64,Float64}[]
    for y = 0:dy:ymax
        R = rstart + im * y; Y = yasym(R); mx = 0.0; at = rstart
        while real(R) > rmin
            k1 = f(R, Y); k2 = f(R - h / 2, Y - h / 2 * k1); k3 = f(R - h / 2, Y - h / 2 * k2); k4 = f(R - h, Y - h * k3)
            Y -= h / 6 * (k1 + 2k2 + 2k3 + k4); R -= h
            nY = norm(Y)
            nY > mx && (mx = nY; at = real(R))
            nY < 1e4 || break
        end
        push!(lines, (y, mx, at))
    end
    lines
end
