# Phase 4: Fock-constrained trial space. Ψ = c0 e^{-εs} f0 + a e^{-εs} D + Σ_{n≥3} terms; a21 = a / c0.
include("../claude3/atoms3.jl")

# seed f0 = e^{εs}[1 + Rψ10 + R²(ψ20|_{a21=0} + ψ21 log R)] truncated at total degree 2, in (s, x, y) terms.
function seed(ε; E = E0)
    Zb = big(Z); πb = BigFloat(π); κb = (πb - 2) / (3πb)
    T1 = Ang((0, 1, :one) => big(1) / 2, (0, 0, :one) => -Zb + ε)             # s^1: Rψ10/s + ε
    q = Dict((0, 0) => big(1) / 2, (2, 0) => big(1) / 2)                         # R²/s² = (1+x²)/2
    Dh = Dict((0, 0) => big(1) / 2, (2, 0) => big(1) / 2, (0, 2) => -big(1))     # D/s²
    T2 = Ang()
    add(k, c) = (T2[k] = get(T2, k, zero(BigFloat)) + c)
    for ((a, b), c) in q
        add((a, b, :one), (1 - 2E) / 12 * c + Zb^2 * c / 3)
        add((a, b, :chi), Zb * c)                                                  # Z q χ
    end
    add((0, 0, :one), Zb^2 / 4); add((2, 0, :one), -Zb^2 / 4)                    # Z² r1 r2 / s² = Z²(1-x²)/4
    add((0, 1, :one), -Zb / 3)                                                   # -Z S r12 / 3 / s²
    for ((a, b), c) in Dh
        add((a, b, :one), -Zb * (πb + 4) / (9πb) * c)                           # harmonic part with a21 = 0
        add((a, b, :lR), -Zb * κb * c)                                            # -Zκ D log h_x
    end
    add((0, 1, :one), ε / 2); add((0, 0, :one), -ε * Zb + ε^2 / 2)              # ε s Rψ10 + ε² s²/2
    T2log = Ang((a, b, :one) => -Zb * κb * c for ((a, b), c) in Dh)              # -Zκ D log s
    [Term(0, 0, mono(0, 0)), Term(1, 0, T1), Term(2, 0, T2), Term(2, 1, T2log)]
end
Dterm() = [Term(2, 0, Ang((0, 0, :one) => big(1) / 2, (2, 0, :one) => big(1) / 2, (0, 2, :one) => -big(1)))]

function build4(N, variant; ε)
    rest = filter(f -> all(t -> t.n >= 3, f), build3(N, variant; ε))
    vcat([seed(ε), Dterm()], rest)
end

# jet check: e^{-εs} f0(s,x,y) vs the repo's 1 + Rψ10 + R²(ψ20 + ψ21 log R) with a21 = 0, along a ray
function jetcheck(ε)
    f = seed(ε)
    x, y = big"0.3", big"0.55"
    r1, r2 = (1 + x) / 2, (1 - x) / 2; R = sqrt(r1^2 + r2^2)
    α = 2atan(r2 / r1); θ = theta(x, y)
    atomv = Dict(:one => big(1), :chi => chi_xy(x, y), :lR => lR_xy(x, y))
    for s in (big"1e-2", big"1e-3", big"1e-4")
        val = sum(s^t.n * log(s)^t.p * sum(c * x^a * y^b * atomv[g] for ((a, b, g), c) in t.Y) for t in f) * exp(-ε * s)
        Rr = s * R
        fock = 1 + Rr * psi10(α, θ; Z = big(2)) + Rr^2 * (psi20(α, θ; Z = big(2), E = E0, a21 = big(0)) + psi21(α, θ; Z = big(2)) * log(Rr))
        @printf("jet check s=%.0e: (seed - Fock jet)/s^3 = %.6f\n", Float64(s), Float64((val - fock) / s^3))
    end
end
