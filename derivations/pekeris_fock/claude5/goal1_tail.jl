# Goal 1: (i) tail sensitivity of the Miller solution (influence of arbitrary data in shells N+1..N+3 on low shells),
# (ii) Miller(E0) vs Ritz eigenvector, (iii) zero of g_N(E) vs Ritz eigenvalue.
include("goal1_backward.jl")
using Random
Random.seed!(7)
open("goal1_tail.out", "w") do io
    for N in (8, 10, 12, 14, 16, 18)
        cols, x0, g0 = miller(E0, N)
        tcols = [(l, m, N1 - l - m) for N1 = N+1:N+3 for l = 0:N1 for m = l:N1-l]
        last = sqrt(sum(weight(t) * abs2(x0[i]) for (i, t) in enumerate(cols) if shell(t) == N))   # size of the last shell
        tv = BigFloat.(randn(length(tcols))) .* last
        _, x1, g1 = miller(E0, N; tail = (tcols, tv))
        dx = shellnorms(cols, x1 - x0, N)
        line = @sprintf("N=%2d tail scale %.2e: |dx| shell 0..4: %s ... shell N: %.2e   dg=%.2e", N, Float64(last),
            join([@sprintf("%.2e", Float64(d)) for d in dx[1:5]], " "), Float64(dx[end]), Float64(abs(g1 - g0)))
        println(line); println(io, line)
        # Ritz vector comparison and eigenvalue of the square truncation vs zero of g_N
        _, v, Er = ritz(N)
        dv = shellnorms(cols, x0 - v, 4)
        # secant solve g_N(E) = 0 near E0
        Ea, Eb = E0 - big"1e-4", E0 + big"1e-4"; ga = miller(Ea, N)[3]; gb = miller(Eb, N)[3]
        for _ = 1:30
            Ec = Eb - gb * (Eb - Ea) / (gb - ga); Ea, ga = Eb, gb; Eb = Ec; gb = miller(Eb, N)[3]
            abs(Eb - Ea) < big"1e-40" && break
        end
        line = @sprintf("      Miller(E0) - Ritz vector, shells 0..4: %s   zero of g_N: %.15f   Ritz E_N: %.15f",
            join([@sprintf("%.2e", Float64(d)) for d in dv], " "), Float64(Eb), Float64(Er))
        println(line); println(io, line); flush(io)
    end
end
