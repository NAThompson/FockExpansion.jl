# Driver: energies (and Gram condition numbers) vs basis size for the variants of atoms_ritz.jl.
include("atoms_ritz.jl")
setprecision(BigFloat, parse(Int, get(ENV, "BITS", "320")))
const ε = big"1.704"
Nmax = parse(Int, get(ARGS, 1, "8"))
variants = Symbol.(split(get(ARGS, 2, "none,log1,log2,atoms2,atoms4"), ","))
hq = BigFloat(1) / parse(Int, get(ENV, "HINV", "16")); kmax = parse(Int, get(ENV, "KMAX", "72"))
t0 = time()
nd = nodes(hq, kmax)
atoms = Dict{Symbol,NTuple{3,Vector{BigFloat}}}()
atoms[:one] = (ones(BigFloat, length(nd.x)), zeros(BigFloat, length(nd.x)), zeros(BigFloat, length(nd.x)))
any(v -> v in (:atoms2, :atoms4), variants) && (atoms[:chi] = atom_table(chi_xy, nd); atoms[:lR] = atom_table(lR_xy, nd))
if :atoms4 in variants
    atoms[:psi30h] = atom_table(psi30h_xy, nd; float64 = true)
    atoms[:lR2] = atom_table((x, y) -> lR_xy(x, y)^2, nd)
end
@printf("nodes %d, atoms %s, %.1fs\\n", length(nd.x), join(keys(atoms), ","), time() - t0)
# quadrature check on a polynomial moment: ∫_{T+} x^2 y^3 = 1/((3)(2+3+2))... exact: ∫0^1 y^3 ∫0^y x^2 = ∫ y^6/3 = 1/21
M = moments(nd, atoms, 2Nmax + 6, 2Nmax + 6)
@printf("moment check: %.3e   (%.1fs)\\n", Float64(M.tab[(:one, :v, :one, :v)][3, 4] / 2 - big(1) / 21), time() - t0)
R = Radial(ε)
open(get(ENV, "OUT", "run_atoms.out"), "a") do io
    for N = 2:Nmax, v in variants
        B = build(N, v; ε)
        t = @elapsed H, S = matrices(B, M, R, ε)
        E, x = ground(H, S)
        c = condest(S)
        line = @sprintf("%-7s N=%2d size=%4d  E-E0=% .4e  cond(S)=%.2e  (%.1fs)", v, N, length(B), Float64(E - E0), Float64(c), t)
        println(line); println(io, line); flush(io)
        # export state, normalized to Ψ(0) = 1: lines "n p a b atom coef"
        isdir("states") || mkdir("states")
        c0 = x[1]  # first basis function carries the constant term 1
        open("states/$(v)_N$(N).txt", "w") do f
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
