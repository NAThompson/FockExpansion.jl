# Export the :fixedlog Rayleigh–Ritz states (normalized Ψ(0) = 1) for independent a21 functionals.
# Format per line: a b c p coefficient   meaning  Ψ = exp(-(u+v+w)/2) Σ coef u^a v^b w^c log(u+v+w)^p.
include("galerkin_log.jl")
isdir("states") || mkdir("states")
for N = 6:13
    r = solve(N, :fixedlog)
    T = Term()
    for (f, c) in zip(r.B, r.x), (k, x) in f
        addto!(T, k, c * x)
    end
    c0 = T[(0, 0, 0, 0, 0)]
    open("states/fixedlog_N$(N).txt", "w") do io
        println(io, "# eps = 1.704  Z = 2  lambda = 1  E_variational = ", r.E)
        for (k, x) in sort(collect(T))
            println(io, k[1], " ", k[2], " ", k[3], " ", k[4], " ", string(x / c0))
        end
    end
end
