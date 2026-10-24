# The hypothesis as stated: one formal solution per selected root,
#   Phi_i = exp(lambda R) R^(-mu) sum_n h_n R^(-n),
#   P h_n + (C - (n-1+mu) Q) h_{n-1} + ((n-2+mu)(n-6+mu) I - Lambda) h_{n-2} = 0,
#   P = lambda^2 + 2 lambda U + T, Q = 2(lambda + U), C = 5 lambda + D,
#   mu = w'Cv / w'Qv (left/right null vectors of P).
# Borel coefficients g_n = h_n / Gamma(n + mu); their radius vs nearest other root.
include("exterior_lib.jl")
using Printf, SpecialFunctions
setprecision(BigFloat, 512)

function formal_solution(b, E, lam, N)
    m = length(b.K)
    U = big.(b.U); Dm = big.(b.D); T = big.(b.W) + 2big(E) * I; Lam = Diagonal(big.(b.lam))
    # refine the root by Newton on det P via smallest singular direction
    l = big(lam)
    for _ = 1:30
        P = l^2 * I + 2l * U + T
        F = svd(Float64.(P))
        v = big.(F.V[:, end]); w = big.(F.U[:, end])
        # inverse iteration in BigFloat for the null vectors
        for _ = 1:4
            Ps = P + big(1e-120) * I; v = Ps \ v; v /= norm(v); w = transpose(Ps) \ w; w /= norm(w)
        end
        dl = -dot(w, P * v) / dot(w, (2l * I + 2U) * v)
        l += dl
        abs(dl) < big(2.0)^-480 && break
    end
    P = l^2 * I + 2l * U + T; Q = 2(l * I + U); C = 5l * I + Dm
    v = big.(randn(m)); w = big.(randn(m))
    for _ = 1:8
        Ps = P + big(1e-120) * I; v = Ps \ v; v /= norm(v); w = transpose(Ps) \ w; w /= norm(w)
    end
    mu = dot(w, C * v) / dot(w, Q * v)
    h = [v]
    # n = 1: P h1 = -(C - mu Q) v, with v'h1 = 0 (bordered)
    Bd = [P v; transpose(v) 0]
    sol = Bd \ vcat(-(C - mu * Q) * v, 0)
    push!(h, sol[1:m])
    for n = 2:N
        A = C - (n - 1 + mu) * Q
        rhs = -A * h[n] - ((n - 2 + mu) * (n - 6 + mu) * I - Lam) * h[n-1]
        Bd = [P A*v; transpose(v) 0]
        sol = Bd \ vcat(rhs, 0)
        h[n] += sol[m+1] * v          # fix the free null component of h_{n-1}
        push!(h, sol[1:m])
    end
    l, mu, h
end

for d in (2, 4, 6)
    b = cuspmodel(d); E = cuspsolve(b).E
    slow, fast = pencil_roots(b, E)
    roots = real(vcat(slow, fast))
    println("d=$d  E=$E")
    for lam in real(slow)
        l, mu, h = formal_solution(b, E, lam, 120)
        g = [norm(h[n+1]) / gamma(n + mu) for n = 0:120]
        rad = Float64(g[101] / g[121])^(1 / 20)   # root-test estimate of the Borel radius
        near = minimum(abs(lam - r) for r in roots if abs(lam - r) > 1e-9)
        @printf("  lambda=%+.6f  mu=%+.6f  Borel radius est=%.4f  nearest other root=%.4f\n", Float64(l), Float64(mu), rad, near)
    end
    flush(stdout)
end
