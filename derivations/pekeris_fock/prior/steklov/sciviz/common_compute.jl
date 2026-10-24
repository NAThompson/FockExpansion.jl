# Shared helpers for the steklov sciviz compute scripts (steklov env, read-only use of the study).
# Run from any directory with
#     julia -t 4 --project=<steklov dir> <script>.jl
const STK = normpath(joinpath(@__DIR__, ".."))
include(joinpath(STK, "production.jl"))
using Serialization

const DATA = joinpath(@__DIR__, "data")
mkpath(DATA)
const TD = Double64

# Cached exterior Steklov forms, read only (get_ext would write missing files into the study).
function load_ext(L, Rs, E, name)
    file = joinpath(STK, "cache", "ext_L$(L)_R$(Rs)_E$(E)_$(name).jls")
    isfile(file) || error("no cached exterior form $file")
    return deserialize(file)
end

# Same assembly as exterior_form_dd (production defaults), but the blocks are kept so that the
# Ritz correction coefficients b = -Abb \ Abt t can be recovered for a given trace t.
function exterior_blocks(model0, E, Rs, nL; zsets, nchan, lchan, nR = 40, nphi = 24, nu = 20,
                         beta = 2.5, edges = [0.25, 0.75, 2.0, 5.0], npanel = 16)
    Z, g = model0.Z, model0.g
    terms = Term[]
    for (za, zb, om) in zsets
        append!(terms, hylleraas_terms(om, za, zb))
    end
    nc = count(harmonic_degree(model0.basis) .<= lchan)
    rule2 = exterior_rule_2d(Rs, TD; nR, nphi, beta, edges, npanel)
    rule3 = exterior_rule_3d(Rs, TD; nR, nphi, nu, beta, edges, npanel)
    Ahh = correction_block(terms, rule2, Rs, E, Z, g, TD)
    A33, A3h = blocks3d_fast(model0.basis, nL, nc, nchan, terms, rule3, Rs, E, Z, g, TD)
    Att = A33[1:nL, 1:nL]
    Atb = hcat(A33[1:nL, (nL+1):end], A3h[1:nL, :])
    Abb = [A33[(nL+1):end, (nL+1):end] A3h[(nL+1):end, :]; A3h[(nL+1):end, :]' Ahh]
    sc = 1 ./ sqrt.(abs.(diag(Abb)))
    ch = cholesky(Symmetric(sc .* Abb .* sc'))
    X = ch.L \ (sc .* Atb')
    S = Att - X' * X
    return (; Att, Atb, Abb, sc, ch, S = (S + S') / 2, terms, nL, nc, nchan,
            basis = model0.basis, Rs)
end

# Ritz corrections for trace coefficients t: b = -Abb^{-1} Abt t (Double64).
ritz_corrections(blk, t) = -(blk.sc .* (blk.ch \ (blk.sc .* (blk.Atb' * TD.(t)))))

# Exterior Ritz solution sum_i t_i exp(S) Y_i + sum_j b_j phi_j at points (r1, r2, r12):
# value and radial derivative d/dR, evaluated in Double64.
function exterior_eval(blk, t, b, pts)
    nL, nc, P = blk.nL, blk.nc, blk.nchan
    n3 = nL + nc * P
    nb = length(blk.terms)
    c3 = vcat(TD.(t), b[1:(nc*P)])
    ch = b[(nc*P+1):end]
    val = zeros(length(pts))
    dR = zeros(length(pts))
    chunks = collect(Iterators.partition(eachindex(pts), cld(length(pts), 8Threads.nthreads())))
    Threads.@threads for chunk in chunks
        V = zeros(TD, 1, n3, 4)
        H = zeros(TD, 1, nb, 4)
        for k in chunk
            x1, x2, x12 = TD.(pts[k])
            eval3d!(V, 1, blk.basis, nL, nc, P, x1, x2, x12, blk.Rs, 2.0, 1.0)
            hylleraas_eval!(H, 1, blk.terms, x1, x2, x12, blk.Rs)
            f = zero(TD)
            g = [zero(TD), zero(TD), zero(TD)]
            for j = 1:n3
                f += c3[j] * V[1, j, 1]
                for c = 1:3
                    g[c] += c3[j] * V[1, j, c+1]
                end
            end
            for j = 1:nb
                f += ch[j] * H[1, j, 1]
                for c = 1:3
                    g[c] += ch[j] * H[1, j, c+1]
                end
            end
            R = sqrt(x1^2 + x2^2)
            val[k] = Float64(f)
            dR[k] = Float64((x1 * g[1] + x2 * g[2] + x12 * g[3]) / R)
        end
    end
    return val, dR
end

# Hypersphere angles (repo convention): r1/R = cos(α/2), r2/R = sin(α/2),
# r12/R = sqrt(1 - sin α cos θ). Steklov coordinates: c = cos(phi), s = sin(phi), u = r12/R.
angles_to_csu(α, θ) = (cos(α / 2), sin(α / 2), sqrt(max(1 - sin(α) * cos(θ), 0.0)))

# Same open grid as common.jl's interior().
open_grid(n; ε = 1e-3) = collect(range(ε, π - ε; length = n))

# Harmonic values Y_i(q, s) for the first `ncol` harmonics at points (c, s, u).
function harmonics_at(basis, csu)
    q = [c^2 - s^2 for (c, s, u) in csu]
    sh = [1 - u^2 for (c, s, u) in csu]
    return harmonic_values(basis, q, sh)
end

U_at(csu) = [cusp_U(c, s, u; Z = 2.0, g = 1.0) for (c, s, u) in csu]
