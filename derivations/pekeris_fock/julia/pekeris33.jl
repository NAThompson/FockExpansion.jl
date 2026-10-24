# Stage 4(a): truncated 33-term Pekeris recurrence (eps^2 = -E), symmetric sector A(l,m,n) = A(m,l,n),
# shells l+m+n <= ω. Solves det(A0 + eps*A1) = 0 for the largest eps (generalized eigenproblem),
# Float64 eigen for a start and BigFloat inverse iteration to refine.
using LinearAlgebra, Printf
include("generated.jl")

function index_set(ω)
    idx = Tuple{Int,Int,Int}[]
    for N = 0:ω, l = 0:N, m = l:N-l
        push!(idx, (l, m, N - l - m))
    end
    idx
end

canon(t) = t[1] <= t[2] ? t : (t[2], t[1], t[3])

function pencil(ω, ::Type{T}; Z = 2, lam = 1) where {T}
    idx = index_set(ω)
    pos = Dict(t => i for (i, t) in enumerate(idx))
    n = length(idx)
    A0 = zeros(T, n, n)
    A1 = zeros(T, n, n)
    for (i, (l, m, k)) in enumerate(idx), (off, f0, f1) in PEKERIS33
        s = canon((l + off[1], m + off[2], k + off[3]))
        (minimum(s) < 0 || !haskey(pos, s)) && continue
        j = pos[s]
        A0[i, j] += T(f0(l, m, k, Z, lam))
        A1[i, j] += T(f1(l, m, k, Z, lam))
    end
    A0, A1
end

# largest real eps with det(A0 + eps A1) = 0
function largest_eps(ω)
    A0, A1 = pencil(ω, Float64)
    λ = eigvals(A0, -A1)
    maximum(real(x) for x in λ if isfinite(x) && abs(imag(x)) < 1e-8)
end

function refine_eps(ω, eps0; bits = 256, iters = 40)
    setprecision(BigFloat, bits) do
        A0, A1 = pencil(ω, BigFloat)
        e = BigFloat(eps0)
        x = ones(BigFloat, size(A0, 1))
        for _ = 1:iters  # Newton on the smallest singular direction via inverse iteration
            F = lu(A0 + e * A1)
            y = F \ (A1 * x)
            x = y / norm(y)
            # Rayleigh-type update: (A0 + e A1) x ≈ -δ A1 x
            r = (A0 + e * A1) * x
            q = A1 * x
            δ = -dot(q, r) / dot(q, q)
            e += δ
            abs(δ) < BigFloat(10)^(-60) && break
        end
        e
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    E0 = big"-2.903724377034119598311"
    println(" ω   size   eps (Float64)          eps (BigFloat)                 E = -eps^2              E - E0")
    for ω = 2:16
        e = largest_eps(ω)
        eb = refine_eps(ω, e)
        E = -eb^2
        @printf("%2d %6d   %.17f   %.22e   %.20e   %.3e\n", ω, length(index_set(ω)), e, eb, E, Float64(E - E0))
    end
end
