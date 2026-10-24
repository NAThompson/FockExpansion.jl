# Production-point fields: matched trace, interior/exterior radial derivatives, mismatch,
# the same for a wrong seed (second matching eigenvector), and Ψ(R, α, θ) on 0.05 ≤ R ≤ 4.
# Output: data/fields.jls (plain arrays).
include(joinpath(@__DIR__, "common_compute.jl"))

const RS, L, D, NAME = 1.0, 10, 32, "B10P6"
const E = E_HE
const NANG_FIELD = 129
const NANG_PSI = 97

timings = Dict{String,Float64}()

# --- Exterior blocks (Double64), checked against the cached Schur complement.
zsets, nchan = basis_spec(NAME)
model0 = channel_model(L)
nL = length(model0.basis.K)
timings["exterior_assembly"] = @elapsed blk = exterior_blocks(model0, E, RS, nL; zsets, nchan, lchan = L)
Sext = Float64.(blk.S)
cached = load_ext(L, RS, E, NAME)
@printf("exterior: nb=%d, |S - S_cached|/|S| = %.2e (%.1fs)\n", size(blk.Abb, 1),
        norm(Sext - Float64.(cached.S)) / norm(Sext), timings["exterior_assembly"])

# --- Interior and matching eigenproblem.
model = channel_model(D)
N = max(44, 2D + 28)
timings["interior_setup"] = @elapsed setup = interior_setup(model, E, RS, L; N)
m = nchannels(model)
eg = eigen(setup.M \ (setup.Bint + Sext))
order = sortperm(abs.(eg.values))
spectrum = eg.values[order]
@printf("lambda_1 = %+.6e, lambda_2 = %+.6e (imag %.1e, %.1e)\n", real(spectrum[1]),
        real(spectrum[2]), imag(spectrum[1]), imag(spectrum[2]))
ref = match(setup, Sext)
pc = physical_constants(setup, ref.a)
@printf("a21 = %.13f (log at L=10 d=32 B10P6: 0.4767478964066)\n", pc.a21)

function seeds_for(v)
    t = real(v)
    a = setup.F \ vcat(t, zeros(m - setup.nL))
    return t, a
end
t1, a1 = seeds_for(eg.vectors[:, order[1]])
scale = a1[1]
t1 ./= scale
a1 ./= scale
# Wrong seed: second-smallest |λ|, scaled to the matched trace's M-norm.
abs(imag(spectrum[2])) < 1e-10 * abs(spectrum[2]) || @warn "second eigenvalue complex"
t2, a2 = seeds_for(eg.vectors[:, order[2]])
s2 = sqrt(dot(t1, setup.M * t1) / dot(t2, setup.M * t2)) * sign(t2[1] == 0 ? 1.0 : t2[1])
t2 .*= s2
a2 .*= s2
@printf("a check: |a1 - match.a| = %.2e\n", norm(a1 - ref.a, Inf))

timings["ritz_corrections"] = @elapsed begin
    b1 = ritz_corrections(blk, t1)
    b2 = ritz_corrections(blk, t2)
end

# --- Fields on the sphere R = R*.
αf = open_grid(NANG_FIELD)
θf = open_grid(NANG_FIELD)
csu_f = [angles_to_csu(a, t) for a in αf, t in θf]      # [iα, iθ]
Yf = harmonics_at(model.basis, vec(csu_f))
Uf = U_at(vec(csu_f))
pts_f = [(RS * c, RS * s, RS * u) for (c, s, u) in vec(csu_f)]

function sphere_fields(t, a, b)
    Phi = Yf * (setup.F * a)              # cusp-gauge trace (all channels; equals Σ t_i Y_i)
    dPhi = Yf * (setup.Fp * a)
    ex = exp.(RS .* Uf)
    Psi = ex .* Phi
    dPsi_int = ex .* (Uf .* Phi .+ dPhi)
    Psi_ext, dPsi_ext = exterior_eval(blk, t, b, pts_f)
    sh = size(csu_f)
    r(x) = reshape(x, sh)
    return (; Phi = r(Phi), Psi = r(Psi), dPsi_int = r(dPsi_int), Psi_ext = r(Psi_ext),
            dPsi_ext = r(dPsi_ext), Phi_trace_L = r(Yf[:, 1:nL] * t))
end
timings["sphere_fields"] = @elapsed begin
    F1 = sphere_fields(t1, a1, b1)
    F2 = sphere_fields(t2, a2, b2)
end
@printf("matched: max|Psi_int - Psi_ext| on sphere = %.2e, max|Phi - Σt Y| = %.2e\n",
        maximum(abs, F1.Psi - F1.Psi_ext), maximum(abs, F1.Phi - F1.Phi_trace_L))

# --- Projections of the fluxes onto the trace space T_L (sphere quadrature).
rule = sphere_rule(64, 64)
csu_q = collect(zip(rule.c, rule.s, rule.u))
Yq = harmonics_at(model.basis, csu_q)
Uq = U_at(csu_q)
gq = exp.(RS .* Uq) .* Yq[:, 1:nL]                 # g_s = exp(R* U) Y_s
pts_q = [(RS * c, RS * s, RS * u) for (c, s, u) in csu_q]
function fluxes(t, a, b)
    Phi = Yq * (setup.F * a)
    dPhi = Yq * (setup.Fp * a)
    dint = exp.(RS .* Uq) .* (Uq .* Phi .+ dPhi)
    _, dext = exterior_eval(blk, t, b, pts_q)
    W = RS^5 .* rule.w
    return (; int_pointwise = gq' * (W .* dint), int_weak = setup.Bint * t,
            ext_pointwise = gq' * (W .* dext), ext_weak = -(Sext * t),
            lambda_M_t = NaN .* t, mism_pointwise = gq' * (W .* (dint .- dext)))
end
P1 = fluxes(t1, a1, b1)
P2 = fluxes(t2, a2, b2)
P1 = merge(P1, (; lambda_M_t = real(spectrum[1]) .* (setup.M * t1)))
P2 = merge(P2, (; lambda_M_t = real(spectrum[2]) .* (setup.M * t2)))
for (nm, P) in (("matched", P1), ("wrong", P2))
    @printf("%s: |int_pw - int_weak| = %.2e, |ext_pw - ext_weak| = %.2e, |mism_pw| = %.2e, |λMt| = %.2e\n",
            nm, norm(P.int_pointwise - P.int_weak, Inf), norm(P.ext_pointwise - P.ext_weak, Inf),
            norm(P.mism_pointwise, Inf), norm(P.lambda_M_t, Inf))
end

# --- Ψ(R, α, θ): interior series for R ≤ R*, exterior Ritz for R > R*.
αp = open_grid(NANG_PSI)
θp = open_grid(NANG_PSI)
csu_p = [angles_to_csu(a, t) for a in αp, t in θp]
Yp = harmonics_at(model.basis, vec(csu_p))
Up = U_at(vec(csu_p))
R_int = collect(0.05:0.05:1.0)
R_ext = collect(1.1:0.1:4.0)
Rall = vcat(R_int, R_ext)
side = vcat(fill(-1, length(R_int)), fill(1, length(R_ext)))
Psi3 = zeros(NANG_PSI, NANG_PSI, length(Rall))
Phi3 = zeros(size(Psi3))
timings["interior_volume"] = @elapsed begin
    ser = fock_series(model, E, N, R_int; nstore = 0)
    for (k, R) in enumerate(R_int)
        Phi = Yp * (ser.F[k] * a1)
        Phi3[:, :, k] = reshape(Phi, NANG_PSI, NANG_PSI)
        Psi3[:, :, k] = reshape(exp.(R .* Up) .* Phi, NANG_PSI, NANG_PSI)
    end
end
timings["exterior_volume"] = @elapsed for (j, R) in enumerate(R_ext)
    k = length(R_int) + j
    pts = [(R * c, R * s, R * u) for (c, s, u) in vec(csu_p)]
    v, _ = exterior_eval(blk, t1, b1, pts)
    Psi3[:, :, k] = reshape(v, NANG_PSI, NANG_PSI)
    Phi3[:, :, k] = reshape(exp.(-R .* Up) .* v, NANG_PSI, NANG_PSI)
end
# Interior series continued a little past R* (where it is still convergent), for comparison.
R_cont = collect(1.0:0.05:1.5)
ser_c = fock_series(model, E, N, R_cont; nstore = 0)
cont_diff = zeros(NANG_PSI, NANG_PSI, length(R_cont))
for (k, R) in enumerate(R_cont)
    pts = [(R * c, R * s, R * u) for (c, s, u) in vec(csu_p)]
    v, _ = exterior_eval(blk, t1, b1, pts)
    Psi_c = exp.(R .* Up) .* (Yp * (ser_c.F[k] * a1))
    cont_diff[:, :, k] = reshape(Psi_c .- v, NANG_PSI, NANG_PSI)
end
println("timings: ", timings)

serialize(joinpath(DATA, "fields.jls"), Dict{String,Any}(
    "Rs" => RS, "L" => L, "d" => D, "basis" => NAME, "E" => E, "N" => N,
    "nL" => nL, "nchannels" => m, "nb_exterior" => size(blk.Abb, 1),
    "labels" => [collect(x) for x in model.basis.labels[1:nL]],
    "spectrum_re" => real.(spectrum), "spectrum_im" => imag.(spectrum),
    "a_matched" => a1[1:min(m, 15)], "a_wrong" => a2[1:min(m, 15)], "a21" => pc.a21,
    "t_matched" => t1, "t_wrong" => t2,
    "alpha" => αf, "theta" => θf,
    "matched" => Dict(String(k) => v for (k, v) in pairs(F1)),
    "wrong" => Dict(String(k) => v for (k, v) in pairs(F2)),
    "proj_matched" => Dict(String(k) => v for (k, v) in pairs(P1)),
    "proj_wrong" => Dict(String(k) => v for (k, v) in pairs(P2)),
    "alpha_psi" => αp, "theta_psi" => θp, "R" => Rall, "side" => side,
    "Psi" => Psi3, "Phi" => Phi3, "R_cont" => R_cont, "cont_diff" => cont_diff,
    "timings" => timings,
))
println("wrote ", joinpath(DATA, "fields.jls"))
