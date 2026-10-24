# Goal 1 numerics: Miller-type backward truncation of the 57-term fixed-eps Pekeris recurrence (singlet sector).
# For energy E (real or complex): unknowns in shells <= N, A(0,0,0) = 1, all equations of shells <= N except the (0,0,0) row.
# x_N(E) -> "minimal solution" if it converges as N grows; g_N(E) = dropped row residual (vanishes at eigenvalues).
using LinearAlgebra, Printf
include("generated57.jl")
setprecision(BigFloat, parse(Int, get(ENV, "BITS", "160")))
const EPS = big"1.704"
const E0 = big"-2.903724377034119598311159245194404446696925310"

idx(N) = [(l, m, N1 - l - m) for N1 = 0:N for l = 0:N1 for m = l:N1-l]
canon(t) = t[1] <= t[2] ? t : (t[2], t[1], t[3])
shell(t) = t[1] + t[2] + t[3]
weight(t) = t[1] == t[2] ? 1 : 2   # multiplicity of (l,m,n) and (m,l,n) in the full coefficient vector

"rows = equations at indices rows; cols = unknown indices; K(E) entries"
function kmat(rows, cols, E::T) where {T}
    pos = Dict(t => i for (i, t) in enumerate(cols))
    K = zeros(T, length(rows), length(cols))
    for (i, (l, m, n)) in enumerate(rows), (off, f0, f1) in PEKERIS57
        s = canon((l + off[1], m + off[2], n + off[3]))
        (minimum(s) < 0 || !haskey(pos, s)) && continue
        K[i, pos[s]] += f0(l, m, n, EPS) + E * f1(l, m, n, EPS)
    end
    K
end

function miller(E::T, N; tail = nothing) where {T}
    cols = idx(N); rows = cols[2:end]
    K = kmat(rows, cols, E)
    rhs = -K[:, 1]
    if tail !== nothing   # prescribed data in shells N+1..N+3 (default for a minimal-solution test: zero)
        tcols, tvals = tail
        Kt = kmat(rows, tcols, E)
        rhs -= Kt * tvals
    end
    x = vcat(one(T), K[:, 2:end] \ rhs)
    g = sum(kmat([(0, 0, 0)], cols, E) .* transpose(x))
    cols, x, g
end

shellnorms(cols, x, Nmax) = [sqrt(sum((weight(t) * abs2(x[i]) for (i, t) in enumerate(cols) if shell(t) == k); init = zero(real(eltype(x))))) for k = 0:Nmax]

function ritz(N)  # eigenvector of the square truncation nearest E0, normalized A000 = 1
    cols = idx(N)
    A0 = kmat(cols, cols, zero(BigFloat)); A1 = kmat(cols, cols, one(BigFloat)) - A0
    F = lu(A0 + E0 * A1); v = ones(BigFloat, length(cols)); E = E0
    for _ = 1:60
        w = F \ (A1 * v); v = w / norm(w)
        En = E0 - one(BigFloat) / (dot(A1 * v, v) == 0 ? 1 : 1) * 0  # placeholder, Rayleigh below
    end
    # Rayleigh-type estimate from A0 v + E A1 v = 0 (least squares in E)
    E = -dot(A1 * v, A0 * v) / dot(A1 * v, A1 * v)
    cols, v / v[1], E
end

if abspath(PROGRAM_FILE) == @__FILE__
    Ns = parse.(Int, split(get(ENV, "NS", "6,8,10,12,14,16,18,20"), ","))
    Es = Any[E0, big"-2.5", Complex{BigFloat}(E0, big"0.05"), Complex{BigFloat}(big"-1.0", big"0.5")]
    open("goal1_backward.out", "w") do io
        for E in Es
            ref = nothing; prev = nothing
            xs = Dict{Int,Any}()
            for N in Ns
                cols, x, g = miller(E, N)
                xs[N] = (cols, x)
                line = @sprintf("E=%-24s N=%2d size=%4d  g_N=% .3e%+.3ei", string(round(ComplexF64(E); digits = 6)), N, length(cols), Float64(real(g)), Float64(imag(g)))
                println(line); println(io, line); flush(io)
            end
            # convergence of the low shells (k <= 4) against the largest N
            colsL, xL = xs[Ns[end]]; posL = Dict(t => i for (i, t) in enumerate(colsL))
            for N in Ns[1:end-1]
                cols, x = xs[N]
                d = sqrt(sum((weight(t) * abs2(x[i] - xL[posL[t]]) for (i, t) in enumerate(cols) if shell(t) <= 4); init = zero(BigFloat)))
                line = @sprintf("   low-shell (k<=4) distance to N=%d solution: N=%2d  %.3e", Ns[end], N, Float64(d))
                println(line); println(io, line)
            end
            sn = shellnorms(colsL, xL, Ns[end])
            println(io, "   shell norms at N=$(Ns[end]): ", join([@sprintf("%.2e", Float64(s)) for s in sn], " "))
            flush(io)
        end
    end
end
