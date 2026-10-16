# Helpers for the truncated Fock expansion figures (fock_truncated_*.jl).
#
#     Ψ_K(R, α, θ) = Σ_{k ≤ K} R^k Σ_p (log R)^p ψ_{k,p}(α, θ)
#
# Coefficients at one angle are stored as c[k+1][p+1] = ψ_{k,p}(α, θ).
isdefined(@__MODULE__, :OUT) || include(joinpath(@__DIR__, "common.jl"))
using Printf

const Z, E, A21 = HE.Z, HE.E, HE.a21

"""Angular coefficients ψ_{k,p}(α, θ) for k ≤ kmax (kmax ≤ 4)."""
function coefficients(α, θ; kmax = 3, step40 = 0.125)
    c = [
        [psi00(α, θ)],
        [psi10(α, θ; Z)],
        [psi20(α, θ; Z, E, a21 = A21), psi21(α, θ; Z)],
        [psi30(α, θ; Z, E, a21 = A21), psi31(α, θ; Z)],
    ]
    if kmax ≥ 4
        push!(
            c,
            [
                psi40(α, θ; Z, E, a21 = A21, step = step40),
                psi41(α, θ; Z, E, a21 = A21),
                psi42(α, θ; Z),
            ],
        )
    end
    c[1:(kmax+1)]
end

# The package's closed forms on the fock_truncated (α, θ) grid are written once, by
# hierarchy_compute.jl, to output/hierarchy_package.nc /fock_grid. psi40 throws (its
# quadrature reaches θ = π) for θ within about 0.01 of π, so θ stops at π - 0.03.
const FOCK_KP = [(0, 0), (1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (4, 0), (4, 1), (4, 2)]
const PACKAGE = joinpath(OUT, "hierarchy_package.nc")

"""(α, θ) grid of the fock_truncated_*.jl data: `n` × `n`."""
fock_grid(n = 64) = (interior(n), collect(range(1e-3, π - 0.03; length = n)))

"""Coefficient vectors `c[iα, iθ]` (as `coefficients(α, θ; kmax = 4)`) and the grid,
read from hierarchy_package.nc /fock_grid."""
function load_fock_grid()
    isfile(PACKAGE) || error("$PACKAGE not found: run hierarchy_compute.jl first")
    NCDataset(PACKAGE) do ds
        haskey(ds.group, "fock_grid") ||
            error("$PACKAGE has no /fock_grid group: rerun hierarchy_compute.jl")
        g = ds.group["fock_grid"]
        αs, θs = g["alpha"][:], g["theta"][:]
        P = Dict(kp => g["psi$(kp[1])$(kp[2])"][:, :] for kp in FOCK_KP)
        c(i, j) = [[P[(k, p)][i, j] for p = 0:(k÷2)] for k = 0:4]
        (; αs, θs, coeffs = [c(i, j) for i in eachindex(αs), j in eachindex(θs)])
    end
end

"""R^k Σ_p (log R)^p ψ_{k,p}: the order-k contribution."""
term(c, k, R) = k + 1 > length(c) ? zero(R) : R^k * evalpoly(log(R), c[k+1])

"""Truncated series Ψ_K."""
psiK(c, K, R) = sum(term(c, k, R) for k = 0:K)

"""Hyperspherical potential V = R(1/r12 - Z/r1 - Z/r2) at angle (α, θ)."""
potential(α, θ) = 1 / sqrt(1 - sin(α) * cos(θ)) - Z * (1 / cos(α / 2) + 1 / sin(α / 2))

"""(H - E)Ψ_K, exact from the recurrence. Every coefficient through order K satisfies its
equation, so only the first two missing orders leave a remainder:
(H-E)Ψ_K = R^{K-1} Σ_p (Vψ_{K,p} - Eψ_{K-1,p}) log^p R - E R^K Σ_p ψ_{K,p} log^p R."""
function residual(c, K, R, V)
    L = log(R)
    cK = c[K+1]
    cKm = K ≥ 1 ? c[K] : [zero(R)]
    n = max(length(cK), length(cKm))
    pad(v) = [v; zeros(n - length(v))]
    s = V .* pad(cK) .- E .* pad(cKm)
    R^(K - 1) * evalpoly(L, s) - E * R^K * evalpoly(L, cK)
end

"""(α, θ) of the electron triangle with side lengths r1, r2 and inter-electron angle θ."""
alpha_of(r1, r2) = 2atan(r2, r1)

const KCOLORS = Makie.wong_colors()[[1, 2, 3, 4, 6]]
