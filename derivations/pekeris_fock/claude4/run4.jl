setprecision(BigFloat, parse(Int, get(ENV, "BITS", "384")))
include("constrained.jl")
PSI30_DIGITS[] = 50
const ε = big"1.704"
jetcheck(ε)
Nmax = parse(Int, get(ARGS, 1, "11")); variants = Symbol.(split(get(ARGS, 2, "min5"), ","))
hq = BigFloat(1) / 16; kmax = 72
nd = nodes(hq, kmax); nn = length(nd.x)
using Serialization
atoms = Dict{Symbol,NTuple{3,Vector{BigFloat}}}(:one => (ones(BigFloat, nn), zeros(BigFloat, nn), zeros(BigFloat, nn)))
fun = Dict(:Lt => Lt_xy, :alp => alp_xy, :V => V_xy, :U => U_xy, :chi => chi_xy, :lR => lR_xy, :psi30h => psi30h_big)
need = Set([:chi, :lR, :Lt, :alp, :V, :U]); :min5psi in variants && push!(need, :psi30h)
for g in need
    cache = "../claude3/atomcache_$(g)_h16_k$(kmax)_b$(precision(BigFloat))_p50.jls"
    atoms[g] = isfile(cache) ? deserialize(cache) : (t = atom_table(fun[g], nd); serialize(cache, t); t)
end
M = moments(nd, atoms, 2Nmax + 6, 2Nmax + 6)
R = Radial(ε)
open(get(ENV, "OUT", "run4.out"), "a") do io
    for N = parse(Int, get(ENV, "NMIN", "4")):Nmax, v in variants
        B = build4(N, v; ε)
        H, S = matrices(B, M, R, ε)
        E, x = ground(H, S)
        c = N <= 9 ? condest(S) : big(NaN)
        line = @sprintf("%-8s N=%2d size=%4d  E-E0=% .4e  a21=%.12f  cond(S)=%.2e", v, N, length(B), Float64(E - E0), Float64(x[2] / x[1]), Float64(c))
        println(line); println(io, line); flush(io); flush(stdout)
        isdir("states") || mkdir("states")
        open("states/c$(v)_N$(N).txt", "w") do f
            println(f, "# eps = 1.704  Z = 2  lambda = 1  E_variational = ", E)
            acc = Dict{Tuple{Int,Int,Int,Int,Symbol},BigFloat}()
            for (bf, cb) in zip(B, x), t in bf, ((a, b, g), cy) in t.Y
                k = (t.n, t.p, a, b, g); acc[k] = get(acc, k, zero(BigFloat)) + cb * cy / x[1]
            end
            for (k, cf) in sort(collect(acc); by = first)
                println(f, k[1], " ", k[2], " ", k[3], " ", k[4], " ", k[5], " ", string(cf))
            end
        end
    end
end
