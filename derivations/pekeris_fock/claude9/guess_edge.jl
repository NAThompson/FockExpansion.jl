# ODE guess on edge data by pivoted QR in BigFloat: report best |R_last|/|R_11| over (r, d).
using LinearAlgebra, Printf
setprecision(BigFloat, 256)
name = ARGS[1]; maxunk = parse(Int, ARGS[2]); rmax = parse(Int, get(ARGS, 3, "12"))
rows = Tuple{BigFloat,Vector{BigFloat}}[]
for l in eachline("edge_series.txt")
    p = split(l)
    p[1] == name || continue
    push!(rows, (parse(BigFloat, p[2]), [parse(BigFloat, x) * factorial(big(i - 1)) for (i, x) in enumerate(p[3:end])]))
end
found = false
for tot = 2:maxunk, r = 1:rmax
    tot % (r + 1) == 0 || continue
    d = tot ÷ (r + 1) - 1
    A = [ (2τ - 1)^j * fd[i+1] for (τ, fd) in rows, i = 0:r for j = 0:d]
    A = reshape(A, length(rows), :)
    A ./= [norm(A[:, k]) for k in 1:size(A, 2)]'
    F = qr(A, ColumnNorm())
    dg = abs.(diag(F.R))
    q = Float64(dg[end] / dg[1]); q2 = Float64(dg[end-1] / dg[1])
    if q < 1e-30 && q2 > 1e3 * q
        @printf("%s: ODE r=%d d=%d (%d unknowns): last %.1e prev %.1e\n", name, r, d, tot, q, q2)
        global found = true
        break
    end
end
found || println(name, ": none up to ", maxunk)
