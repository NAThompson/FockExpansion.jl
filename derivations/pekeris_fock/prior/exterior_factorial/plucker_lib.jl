# Pole-free homogeneous (Plucker) representation of the slow subspace.
#
# First-order system X = (Phi, Phi'), X' = A(R) X with
#   A(R) = A0 + A1/R + A2/R^2,
#   A0 = [0 I; -T -2U],  A1 = [0 0; -D -5I],  A2 = [0 0; Lambda 0].
# The slow frame X (2m x m) gives Xi = X_1 ^ ... ^ X_m, which solves the LINEAR
# compound system Xi' = A^[m](R) Xi (singular only at R = 0).
# Formal solution Xi = exp(sigma R) R^beta sum_n xi_n R^{-n}:
#   (C0 - sigma) xi_n = ((beta - n + 1) I - C1) xi_{n-1} - C2 xi_{n-2},
#   beta = w'C1 xi_0 / w'xi_0  (w: left null vector of C0 - sigma).
# Then xi(R) = xi_0 + sum_k B_k/(x)_{k+1}, x = R/omega, as for Y.
include("exterior_lib.jl")

"""All k-subsets of 1:n in lexicographic order."""
function subsets(n, k)
    out = Vector{Vector{Int}}()
    rec(start, cur) =
        length(cur) == k ? push!(out, copy(cur)) :
        for i = start:n
            push!(cur, i); rec(i + 1, cur); pop!(cur)
        end
    rec(1, Int[])
    out
end

"""Plucker vector of an n x k frame: all k x k minors over row subsets S."""
wedge(M, S) = [det(M[s, :]) for s in S]

"""Additive compound A^[k] (linear in A): column T = sum_j wedge(E_T with column j -> A e_{t_j})."""
function additive_compound(A, S)
    n = size(A, 1)
    idx = Dict(s => i for (i, s) in enumerate(S))
    C = zeros(eltype(A), length(S), length(S))
    for (c, t) in enumerate(S)
        E = Matrix{eltype(A)}(I, n, n)[:, t]
        for j in eachindex(t)
            M = copy(E); M[:, j] = A[:, t[j]]
            C[:, c] += wedge(M, S)
        end
    end
    C
end

function first_order_blocks(b, E)
    m = length(b.K)
    U = big.(b.U); Dm = big.(b.D); T = big.(b.W) + 2big(E) * I; Lam = Diagonal(big.(b.lam))
    Z = zeros(BigFloat, m, m); Id = Matrix{BigFloat}(I, m, m)
    A0 = [Z Id; -T -2U]
    A1 = [Z Z; -Dm -5Id]
    A2 = [Z Z; Matrix(Lam) Z]
    A0, A1, A2
end

"""Formal Plucker series xi_0..xi_N (all finalized), sigma, beta, and the A0^[m] spectrum gap."""
function plucker_series(b, E, N; prec = 1024)
    setprecision(BigFloat, prec) do
        m = length(b.K)
        S = subsets(2m, m)
        A0, A1, A2 = first_order_blocks(b, E)
        C0, C1, C2 = (additive_compound(A, S) for A in (A0, A1, A2))
        # slow frame of A0: invariant subspace [I; Y0]
        Y0 = slow_solvent(b, E; prec)
        xi0 = wedge([Matrix{BigFloat}(I, m, m); Y0], S)
        sigma = tr(Y0)
        # left null vector of C0 - sigma by inverse iteration
        Ms = C0 - sigma * I
        w = big.(randn(length(S)))
        F = lu(transpose(Ms) + big(1e-200) * I)
        for _ = 1:6
            w = F \ w; w /= norm(w)
        end
        wx = dot(w, xi0)
        beta = dot(w, C1 * xi0) / wx
        M = lu(Ms + xi0 * transpose(w))   # nonsingular for a simple eigenvalue
        xi = [xi0, zeros(BigFloat, length(S))]
        r = (beta * I - C1) * xi0
        xi[2] = M \ r
        for n = 2:N+1
            r = ((beta - n + 1) * I - C1) * xi[n] - C2 * xi[n-1]
            alpha = dot(w, r) / ((n - 1) * wx)
            xi[n] += alpha * xi0                     # finalize xi_{n-1}
            r += alpha * ((beta - n + 1) * I - C1) * xi0
            push!(xi, M \ r)
        end
        (; xi = xi[1:N+1], sigma, beta, C0, C1, C2, S, resid0 = norm(C0 * xi0 - sigma * xi0))
    end
end

"""Projective sine distance between two vectors."""
function proj_dist(a, b)
    a = a / norm(a); b = b / norm(b)
    norm(a - dot(b, a) * b)
end

"""Recover the graph Y from Xi: Y_ij = s_ij Xi[rows (1:m minus j) + (m+i)] / Xi[1:m]."""
function graph_from_plucker(Xi, S, m)
    idx = Dict(s => k for (k, s) in enumerate(S))
    top = Xi[idx[collect(1:m)]]
    Y = zeros(eltype(Xi), m, m)
    for i = 1:m, j = 1:m
        rows = sort(vcat(setdiff(1:m, j), m + i))
        # [I;Y] minor with row j replaced by row m+i: moving it to position j gives (-1)^(m-j)
        Y[i, j] = (-1)^(m - j) * Xi[idx[rows]] / top
    end
    Y
end

"""Pairing det[P | Q] = sum_S sign(S) det P[S,:] det Q[S^c,:] for 2m x m frames."""
function plucker_pairing(wP, XiQ, S, n)
    idx = Dict(s => k for (k, s) in enumerate(S))
    m = length(S[1])
    acc = zero(promote_type(eltype(wP), eltype(XiQ)))
    for (k, s) in enumerate(S)
        sc = setdiff(1:n, s)
        acc += (-1)^(sum(s) + m * (m + 1) ÷ 2) * wP[k] * XiQ[idx[sc]]
    end
    acc
end

"""factorial_coeffs for vector-valued series (same formula as exterior_lib's matrix version)."""
function factorial_coeffs_vec(xs, omega, K; prec = precision(xs[1][1]), stir = stirling1(K))
    setprecision(BigFloat, prec) do
        w = big(omega)
        [sum((stir[k+1, j+1] * w^(-j - 1)) * xs[j+2] for j = 0:k) for k = 0:K]
    end
end
