# Goal 2 admissibility: local expansion of the psi30 atom f = h_x^3 psi30 near faces and corners of 0 <= x <= y <= 1.
# Fit f(path(tau)) = sum_k a_k tau^k + b_k tau^k log tau + c_k tau^(k+1/2) on tau = 2^-j and report the leading non-analytic term.
setprecision(BigFloat, 320)
include("../claude3/atoms3.jl")
PSI30_DIGITS[] = 70
using Printf, LinearAlgebra
J0 = parse(Int, get(ENV, "J0", "6"))
f(x, y) = chi_xy(x, y)
paths = [("w-face y=1 at x=0.3", τ -> (big"0.3", 1 - τ)),
         ("u-face y=x at x=0.4", τ -> (big"0.4", big"0.4" + τ)),
         ("e-e corner along x=0.5y", τ -> (τ / 2, τ)),
         ("nucleus corner along (1-2t, 1-t)", τ -> (1 - 2τ, 1 - τ))]
K = parse(Int, get(ENV, "K", "6"))
for (name, path) in paths
    τs = [big(2.0)^(-j) for j = J0:J0+3K+6]
    vals = [f(path(τ)...) for τ in τs]
    cols = Function[]
    labels = String[]
    for k = 0:K
        push!(cols, τ -> τ^k); push!(labels, "τ^$k")
        k >= 1 && (push!(cols, τ -> τ^k * log(τ)); push!(labels, "τ^$k log τ"))
        push!(cols, τ -> τ^k * sqrt(τ)); push!(labels, "τ^$(k).5")
    end
    M = [c(τ) for τ in τs, c in cols]
    coef = M \ vals
    res = norm(M * coef - vals)
    println(name, @sprintf("  (fit residual %.1e)", Float64(res)))
    for (l, c) in zip(labels, coef)
        (occursin("log", l) || occursin(".5", l)) && abs(c) > 1e-12 && println("    ", l, " : ", @sprintf("%.6e", Float64(c)))
    end
    flush(stdout)
end
