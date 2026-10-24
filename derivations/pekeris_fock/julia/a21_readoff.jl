# Stage 4(c): crude a21 read-off from the log-augmented Rayleigh–Ritz state (variant :fixedlog).
# psi20(Ω) = [P2 - σ P1/2 + σ^2 P0/8](Ω) + log(σ/R) Lp2(Ω) at R = 1, divided by Ψ(0) = P0;
# a21 = <psi20, v> / <v, v> on S^5 (S-state measure sin^2 α sin θ dα dθ), v = sin α cos θ.
using LinearAlgebra, Printf
using FockExpansion
include("galerkin_log.jl")

function gausslegendre(n)
    J = SymTridiagonal(zeros(n), [k / sqrt(4k^2 - 1) for k = 1:n-1])
    x, V = eigen(J)
    x, 2 .* V[1, :] .^ 2
end

function project_v(f; n = 160)
    x, wq = gausslegendre(n)
    num = 0.0; den = 0.0
    for (xa, wa) in zip(x, wq), (xt, wt) in zip(x, wq)
        α = π * (xa + 1) / 2; θ = π * (xt + 1) / 2
        μ = sin(α)^2 * sin(θ) * wa * wt
        vv = sin(α) * cos(θ)
        num += μ * f(α, θ) * vv
        den += μ * vv^2
    end
    num / den
end

evalterm(t::Term, u, v, w; deg) = sum((Float64(c) * u^k[1] * v^k[2] * w^k[3] for (k, c) in t if k[4] == 0 && k[1] + k[2] + k[3] == deg); init = 0.0)

function a21_from_state(r, epsv; Z = 2)
    κ = (π - 2) / (3π)
    # collect Taylor pieces of the polynomial part (log-free terms) and the degree-2 log part
    P = Term(); L = Term()
    for (f, c) in zip(r.B, r.x), (k, x) in f
        addto!(k[4] == 0 ? P : L, (k[1], k[2], k[3], 0, 0), c * x)
    end
    c0 = Float64(get(P, (0, 0, 0, 0, 0), big(0)))
    function ψ20(α, θ)
        r1, r2 = cos(α / 2), sin(α / 2)
        r12 = sqrt(max(r1^2 + r2^2 - 2r1 * r2 * cos(θ), 0.0))
        u = epsv * (r2 + r12 - r1); v = epsv * (r1 + r12 - r2); w = 2epsv * (r1 + r2 - r12)
        σ = u + v + w
        p0 = evalterm(P, u, v, w; deg = 0); p1 = evalterm(P, u, v, w; deg = 1); p2 = evalterm(P, u, v, w; deg = 2)
        l2 = evalterm(L, u, v, w; deg = 2)
        (p2 - σ * p1 / 2 + σ^2 * p0 / 8 + log(σ) * l2) / c0   # R = 1 so log(σ/R) = log σ
    end
    project_v(ψ20)
end

if abspath(PROGRAM_FILE) == @__FILE__
    E0 = -2.903724377034119598311
    # convention check: the repo's psi20 with a21 = 0.47674787900 must project back to a21
    @printf("convention check: P_v[psi20(a21=0.476747879)] = %.12f\n",
        project_v((α, θ) -> psi20(α, θ; Z = 2.0, E = E0, a21 = 0.47674787900)))
    epsv = 1.704
    for N = 4:parse(Int, get(ARGS, 1, "12"))
        r = solve(N, :fixedlog)
        @printf("N=%2d size=%4d  E-E0=%.3e  a21 ≈ %.8f\n", N, r.n, Float64(r.E) - E0, a21_from_state(r, epsv))
    end
end
