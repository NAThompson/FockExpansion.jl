# Borel-plane singularity map via Pade poles of a scalar projection of the Borel series.
include("plucker_lib.jl")
using Printf, Random
setprecision(BigFloat, 1024)
b = cuspmodel(2); E = cuspsolve(b).E
N = 130
ps = plucker_series(b, E, N; prec = 1024)
Ys = riccati_coeffs(b, E, N; prec = 1024)

"""Pade [L/L] denominator roots for f(t) = sum_j c_j t^j."""
function pade_poles(c, L)
    # denominator q (q0 = 1): sum_{i=0}^L q_i c_{L+k-i} = 0, k = 1..L
    A = [c[L+k-i+1] for k = 1:L, i = 1:L]
    rhs = [-c[L+k+1] for k = 1:L]
    q = vcat(big(1), A \ rhs)
    aberth(q)
end

"""All roots of sum q_i t^i by Aberth iteration in Complex{BigFloat}."""
function aberth(q; iters = 500)
    n = length(q) - 1
    qc = Complex{BigFloat}.(q)
    pv(z) = foldr((a, acc) -> acc * z + a, qc; init = zero(z))
    dq = [i * qc[i+1] for i = 1:n]
    dv(z) = foldr((a, acc) -> acc * z + a, dq; init = zero(z))
    rad = Float64(abs(q[1] / q[end]))^(1 / n)
    z = [Complex{BigFloat}(rad * cis(2pi * (k + 0.25) / n)) for k = 0:n-1]
    for _ = 1:iters
        maxstep = 0.0
        for k = 1:n
            ratio = pv(z[k]) / dv(z[k])
            s = sum(1 / (z[k] - z[j]) for j = 1:n if j != k)
            w = ratio / (1 - ratio * s)
            z[k] -= w
            maxstep = max(maxstep, Float64(abs(w)))
        end
        maxstep < 1e-40 && break
    end
    ComplexF64.(z)
end

Random.seed!(3)
ell = big.(randn(length(ps.xi[1])))
cx = [dot(ell, ps.xi[j+2]) / factorial(big(j)) for j = 0:N-2]
cy = [Ys[j+2][1, 2] + Ys[j+2][3, 4] for j = 0:N-2] ./ [factorial(big(j)) for j = 0:N-2]
ev = eigvals(Float64.(ps.C0)); cand = sort(real(Float64(ps.sigma) .- ev), by = abs)[2:end]
println("candidate Borel singularities sigma - (other 4-root sums): ", round.(cand[1:12], digits = 3))
open(joinpath(@__DIR__, "plucker_pade.csv"), "w") do io
    println(io, "series,L,re,im")
    for (name, c) in (("plucker", cx), ("Ychart", cy)), L in (40, 60)
        p = pade_poles(c, L)
        foreach(z -> println(io, "$name,$L,$(real(z)),$(imag(z))"), p)
        near = sort(p, by = abs)[1:8]
        println("$name L=$L nearest Pade poles: ", [round(z, digits = 3) for z in near])
    end
end
