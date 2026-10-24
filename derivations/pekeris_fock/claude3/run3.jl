# Phase 3 driver: Ritz energies with Codex's minimal order-2 atoms and a 50-digit psi30 atom; exports states.
setprecision(BigFloat, parse(Int, get(ENV, "BITS", "384")))   # set BEFORE including constants (Codex review)
include("atoms3.jl")
PSI30_DIGITS[] = parse(Int, get(ENV, "PSIDIG", "50"))
const ε = big"1.704"
Nmax = parse(Int, get(ARGS, 1, "11"))
variants = Symbol.(split(get(ARGS, 2, "min5"), ","))
hq = BigFloat(1) / parse(Int, get(ENV, "HINV", "16")); kmax = parse(Int, get(ENV, "KMAX", "72"))
t0 = time()
nd = nodes(hq, kmax)
nn = length(nd.x)
atoms = Dict{Symbol,NTuple{3,Vector{BigFloat}}}(:one => (ones(BigFloat, nn), zeros(BigFloat, nn), zeros(BigFloat, nn)))
need = Set{Symbol}()
for v in variants
    v in (:min5, :min5psi) && union!(need, (:Lt, :alp, :V, :U))
    v in (:atoms2psi,) && union!(need, (:chi, :lR))
    v in (:min5psi, :atoms2psi) && push!(need, :psi30h)
end
fun = Dict(:Lt => Lt_xy, :alp => alp_xy, :V => V_xy, :U => U_xy, :chi => chi_xy, :lR => lR_xy, :psi30h => psi30h_big)
using Serialization
for g in need
    cache = "atomcache_$(g)_h$(get(ENV, "HINV", "16"))_k$(kmax)_b$(precision(BigFloat))_p$(PSI30_DIGITS[]).jls"
    t = @elapsed atoms[g] = isfile(cache) ? deserialize(cache) : (tab = atom_table(fun[g], nd); serialize(cache, tab); tab)
    @printf("atom %s tabulated at %d nodes (%.0fs)\n", g, nn, t); flush(stdout)
end
M = moments(nd, atoms, 2Nmax + 6, 2Nmax + 6)
@printf("moments done (%.0fs)\n", time() - t0); flush(stdout)
R = Radial(ε)
tag = get(ENV, "TAG", "")
open(get(ENV, "OUT", "run3.out"), "a") do io
    for N = parse(Int, get(ENV, "NMIN", "4")):Nmax, v in variants
        B = build3(N, v; ε)
        t = @elapsed H, S = matrices(B, M, R, ε)
        E, x = ground(H, S)
        c = N <= parse(Int, get(ENV, "CONDMAX", "9")) ? condest(S) : big(NaN)  # condest is memory-hungry at large N
        res = norm(H * x - E * (S * x)) / norm(H * x)
        line = @sprintf("%-9s N=%2d size=%4d  E-E0=% .4e  cond(S)=%.2e  resid=%.1e  (%.1fs)", v, N, length(B), Float64(E - E0), Float64(c), Float64(res), t)
        println(line); println(io, line); flush(io); flush(stdout)
        isdir("states") || mkdir("states")
        c0 = x[1]
        open("states/$(v)$(tag)_N$(N).txt", "w") do f
            println(f, "# eps = 1.704  Z = 2  lambda = 1  E_variational = ", E)
            acc = Dict{Tuple{Int,Int,Int,Int,Symbol},BigFloat}()
            for (bf, cb) in zip(B, x), t in bf, ((a, b, g), cy) in t.Y
                k = (t.n, t.p, a, b, g)
                acc[k] = get(acc, k, zero(BigFloat)) + cb * cy / c0
            end
            for (k, cf) in sort(collect(acc); by = first)
                println(f, k[1], " ", k[2], " ", k[3], " ", k[4], " ", k[5], " ", string(cf))
            end
        end
    end
end
