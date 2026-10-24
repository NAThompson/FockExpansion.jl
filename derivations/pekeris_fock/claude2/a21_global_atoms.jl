# Global finite-part a21 estimator (Codex's codex/codex_a21_functional.jl, generalized by Claude to the (s, x, y) atom
# states written by run_atoms.jl: lines "n p a b atom coef", Ψ = e^{-εs} Σ coef s^n log^p s x^a y^b atom(x, y), Ψ(0) = 1).
#   a21 = -D/4 - b/8,  D = FP ∫ P_v[C Ψ / R² - E Ψ / R] dR  (a21_global/FINITE_PART.md; Codex re-derivation).
using LinearAlgebra, Printf
include("atoms_ritz.jl")

function gauss01(n)
    x, V = eigen(SymTridiagonal(zeros(n), [k / sqrt(4k * k - 1) for k = 1:n-1]))
    ((x .+ 1) ./ 2, V[1, :] .^ 2)
end

# ∫_0^∞ s^(h-1) e^{-εs} log^p s ds
function L0(h, p, ε)
    base = factorial(big(h - 1)) / ε^h
    p == 0 && return base
    p0 = polygamma_int(0, h) - log(ε)
    p == 1 && return base * p0
    p1 = polygamma_int(1, h)
    p == 2 && return base * (p0^2 + p1)
    error("p > 2")
end

function read_state(path)
    lines = readlines(path)
    E = parse(BigFloat, split(lines[1], "E_variational = ")[2])
    terms = [(parse.(Int, split(l)[1:4])..., Symbol(split(l)[5]), parse(BigFloat, split(l)[6])) for l in lines[2:end]]
    E, terms
end

const ATOMFUN = Dict(:one => (x, y) -> one(x), :chi => chi_xy, :lR => lR_xy, :lR2 => (x, y) -> lR_xy(x, y)^2,
    :psi30h => (x, y) -> psi30h_big(x, y))

function a21_global(path, nq; bits = 192)
    setprecision(BigFloat, bits) do
        E, terms = read_state(path)
        ε = big"1.704"
        γ = BigFloat(Base.MathConstants.eulergamma)
        xq, wq = gauss01(nq)
        D = zero(BigFloat); normcheck = zero(BigFloat)
        need = unique(t[5] for t in terms)
        for (iy, yf) in enumerate(xq), (it, tf) in enumerate(xq)
            y = BigFloat(yf); x = y * BigFloat(tf); hx = sqrt((1 + x^2) / 2)
            vv = 1 - 2y^2 / (1 + x^2); Vs = -8 / (1 - x^2) + 1 / y
            weight = 64 / BigFloat(π) * y * (1 - x^2) / (1 + x^2)^3 * vv * y * BigFloat(wq[iy]) * BigFloat(wq[it])
            normcheck += weight * vv
            av = Dict(g => ATOMFUN[g](x, y) for g in need)
            dp = zero(BigFloat)
            for (n, p, a, b, g, c) in terms
                T = c * x^a * y^b * av[g]
                if n == 0
                    @assert p == 0
                    dp += T * (Vs * ε * (log(ε) + γ - 1 - log(hx)) + E * (γ + log(ε) - log(hx)))
                elseif n == 1
                    @assert p == 0
                    dp += T * (Vs * (-γ - log(ε) + log(hx)) - E / ε)
                else
                    dp += T * (Vs * L0(n - 1, p, ε) - E * L0(n, p, ε))
                end
            end
            D += weight * dp
        end
        bcoef = -2 * (BigFloat(π) - 2) / (3BigFloat(π))
        -D / 4 - bcoef / 8, normcheck
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    for path in ARGS
        for nq in (32, 48)
            a, nc = a21_global(path, nq)
            @printf("%s nq=%d a21_global=%.12f  P_v[v]=%.12f\n", basename(path), nq, Float64(a), Float64(nc))
            flush(stdout)
        end
    end
end
