# steklov_fields.nc and steklov_wavefunction.nc from data/fields.jls.
#     julia -t auto --project=examples/sciviz <this file>
include("/Users/nathompson7/FockExpansion.jl/examples/sciviz/common.jl")
using Printf, Serialization

const DATA = joinpath(@__DIR__, "data")
F = deserialize(joinpath(DATA, "fields.jls"))
Rs, L, d, basis = F["Rs"], F["L"], F["d"], F["basis"]
setting = "R* = $Rs, L = $L, d = $d, exterior $basis, E = E0"
M, W = F["matched"], F["wrong"]
λ1, λ2 = F["spectrum_re"][1], F["spectrum_re"][2]

# ---------------------------------------------------------------------------------------------
# steklov_fields.nc: angular maps on the matching sphere R = R*
# ---------------------------------------------------------------------------------------------
heat(cm = "auto"; kw...) = view_hint(vis = "heatmap", colormap = cm; kw...)
sym = view_hint(vis = "heatmap", scale = "symlog")
fields = [
    "Phi_trace_cusp" => M["Phi"],
    "Psi_trace" => M["Psi"],
    "dRPsi_interior" => M["dPsi_int"],
    "dRPsi_exterior" => M["dPsi_ext"],
    "dRPsi_int_minus_ext" => M["dPsi_int"] .- M["dPsi_ext"],
    "Psi_int_minus_ext" => M["Psi"] .- M["Psi_ext"],
    "Psi_trace_wrong_seed" => W["Psi"],
    "dRPsi_interior_wrong_seed" => W["dPsi_int"],
    "dRPsi_exterior_wrong_seed" => W["dPsi_ext"],
    "dRPsi_int_minus_ext_wrong_seed" => W["dPsi_int"] .- W["dPsi_ext"],
]
# Mismatch relative to the size of the flux, so both seeds share one colour scale.
relmis(X) = (X["dPsi_int"] .- X["dPsi_ext"]) ./ maximum(abs, X["dPsi_int"])
push!(fields, "dRPsi_int_minus_ext_relative" => relmis(M))
push!(fields, "dRPsi_int_minus_ext_relative_wrong_seed" => relmis(W))

attribs = Dict{String,Any}(
    "Phi_trace_cusp" => ["long_name" => "matched trace, cusp gauge Φ(R*, Ω) = Σ_{deg ≤ L} t_i Y_i", heat()...],
    "Psi_trace" => ["long_name" => "matched trace, physical Ψ(R*, Ω) = exp(R* U) Φ (Ψ(0) = 1)", heat()...],
    "dRPsi_interior" => ["long_name" => "∂_R Ψ at R* from the interior Fock series (all channels)", heat()...],
    "dRPsi_exterior" => ["long_name" => "∂_R Ψ at R* from the exterior Ritz solution (pointwise)", heat()...],
    "dRPsi_int_minus_ext" => ["long_name" => "flux mismatch ∂_R Ψ_int − ∂_R Ψ_ext at R*, matched seeds", sym...],
    "Psi_int_minus_ext" => ["long_name" => "value jump Ψ_int − Ψ_ext at R* (zero by construction: shared trace)", sym...],
    "Psi_trace_wrong_seed" => ["long_name" => "trace for the 2nd matching eigenvector (λ = $(@sprintf("%.3g", λ2))), same M-norm", heat()...],
    "dRPsi_interior_wrong_seed" => ["long_name" => "∂_R Ψ_int at R*, wrong seed (2nd eigenvector)", heat()...],
    "dRPsi_exterior_wrong_seed" => ["long_name" => "∂_R Ψ_ext at R*, wrong seed (2nd eigenvector)", heat()...],
    "dRPsi_int_minus_ext_wrong_seed" => ["long_name" => "flux mismatch at R*, wrong seed: the boundary condition is NOT satisfied", sym...],
    "dRPsi_int_minus_ext_relative" => ["long_name" => "flux mismatch / max|∂_R Ψ_int|, matched", sym...,
                                       view_hint(domain = [-1.0, 1.0])...],
    "dRPsi_int_minus_ext_relative_wrong_seed" => ["long_name" => "flux mismatch / max|∂_R Ψ_int|, wrong seed", sym...,
                                                  view_hint(domain = [-1.0, 1.0])...],
)
path = joinpath(OUT, "steklov_fields.nc")
write_angular(path, F["alpha"], F["theta"], fields; attribs,
              overlaid = ["dRPsi_int_minus_ext", "dRPsi_int_minus_ext_wrong_seed",
                          "dRPsi_int_minus_ext_relative", "dRPsi_int_minus_ext_relative_wrong_seed"],
              global_attribs = [
                  "title" => "Steklov–Poincaré matching for helium on the sphere R = R*",
                  "setting" => setting,
                  "lambda_min" => λ1, "lambda_2" => λ2, "a21" => F["a21"],
                  "nchannels_interior" => F["nchannels"], "n_trace" => F["nL"],
                  "nb_exterior" => F["nb_exterior"],
                  "comment" => "Fields are on an open (α, θ) grid; r1/R = cos(α/2), r2/R = sin(α/2), r12/R = sqrt(1 − sin α cos θ). Ψ normalized by Ψ(0) = 1 (a_1 = 1).",
              ])

# Flux projections onto the trace space T_L and the matching spectrum.
NCDataset(path, "a") do ds
    g = defGroup(ds, "trace_projection")
    nL = F["nL"]
    defVar(g, "harmonic", collect(1:nL), ("harmonic",); attrib = [
        "long_name" => "trace harmonic index i (Y_i, label (m, j), K = 2(m+j))",
        "H5WEB_slice_labels" => "/trace_projection/harmonic_label"])
    defVar(g, "harmonic_label", ["Y($(l[1]),$(l[2]))" for l in F["labels"]], ("harmonic",))
    for (tag, P, λ) in (("matched", F["proj_matched"], λ1), ("wrong_seed", F["proj_wrong"], λ2))
        names = ["int_weak", "ext_weak", "int_pointwise", "ext_pointwise", "mism_pointwise", "lambda_M_t"]
        desc = Dict(
            "int_weak" => "(B_int t)_i: interior flux used by the matching",
            "ext_weak" => "−(A_ext t)_i: exterior flux implied by the Ritz energy form",
            "int_pointwise" => "R*⁵⟨g_i, ∂_R Ψ_int⟩ by sphere quadrature (should equal int_weak)",
            "ext_pointwise" => "R*⁵⟨g_i, ∂_R Ψ_ext⟩ from the pointwise Ritz derivative",
            "mism_pointwise" => "R*⁵⟨g_i, ∂_R Ψ_int − ∂_R Ψ_ext⟩ (pointwise mismatch seen by T_L)",
            "lambda_M_t" => "λ (M t)_i: what the eigen-equation says the projected mismatch is",
        )
        for n in names
            defVar(g, "$(n)_$tag", P[n], ("harmonic",); attrib = [
                "long_name" => "$(desc[n]), $tag", "H5WEB_label" => n])
        end
        defVar(g, "abs_flux_$tag", hcat([abs.(P[n]) for n in names]...), ("harmonic", "quantity_$tag");
               attrib = ["long_name" => "|flux projections| per trace harmonic, $tag (λ = $(@sprintf("%.3e", λ)))",
                         view_hint(vis = "line", scale = "log", x_dim = "harmonic", multi_dim = "quantity_$tag")...])
        defVar(g, "quantity_$tag", collect(1:length(names)), ("quantity_$tag",); attrib = [
            "long_name" => join(names, ", "), "H5WEB_slice_labels" => "/trace_projection/quantity_label_$tag"])
        defVar(g, "quantity_label_$tag", names, ("quantity_$tag",))
    end
    g = defGroup(ds, "spectrum")
    n = length(F["spectrum_re"])
    defVar(g, "index", collect(1:n), ("index",); attrib = ["long_name" => "eigenvalue rank by |λ|"])
    defVar(g, "abs_lambda", abs.(F["spectrum_re"] .+ im .* F["spectrum_im"]), ("index",); attrib = [
        "long_name" => "|λ| of M⁻¹(B_int + A_ext), sorted; λ_1 ≈ 0 is the matched solution",
        view_hint(vis = "line", scale = "log")...])
    defVar(g, "lambda_real", F["spectrum_re"], ("index",); attrib = ["long_name" => "Re λ"])
    defVar(g, "lambda_imag", F["spectrum_im"], ("index",); attrib = ["long_name" => "Im λ"])
    defVar(g, "seeds_matched", F["a_matched"], ("seed",); attrib = [
        "long_name" => "cusp-gauge seeds a_1..a_15 (a_1 = 1), matched"])
    defVar(g, "seeds_wrong_seed", F["a_wrong"], ("seed",); attrib = [
        "long_name" => "seeds a = F⁻¹[t; 0] for the 2nd eigenvector (same trace norm)"])
end
println("wrote ", path)

# ---------------------------------------------------------------------------------------------
# steklov_wavefunction.nc: Ψ(R, α, θ), interior series for R ≤ R*, exterior Ritz beyond
# ---------------------------------------------------------------------------------------------
αp, θp, R = F["alpha_psi"], F["theta_psi"], F["R"]
path = joinpath(OUT, "steklov_wavefunction.nc")
isfile(path) && rm(path)
NCDataset(path, "c"; attrib = [
    "title" => "Glued Steklov wavefunction Ψ(R, α, θ) for helium at E0",
    "setting" => setting,
    "comment" => "R ≤ R*: interior cusp-gauge Fock channel series; R > R*: exterior Double64 Ritz solution. Ψ(0) = 1.",
]) do ds
    defVar(ds, "alpha", αp, ("alpha",); attrib = ["units" => "rad", "long_name" => "α (hyperangle)", PI_TICKS])
    defVar(ds, "theta", θp, ("theta",); attrib = ["units" => "rad", "long_name" => "θ (angle between r₁, r₂)", PI_TICKS])
    defVar(ds, "R", R, ("R",); attrib = ["units" => "bohr", "long_name" => "hyperradius R"])
    xyz = write_xyz!(ds, αp, θp, ("alpha", "theta"))
    angular = write_coalescence!(ds)
    # Julia (alpha, theta, R) is HDF5 (R, theta, alpha): R is the slider.
    slider = view_hint(probe_dim = "R", playback_dim = "R", playback_range = "frame")
    for (name, A, ln) in (
        ("Psi", F["Psi"], "physical Ψ(R, α, θ), slider over R"),
        ("Phi", F["Phi"], "cusp-gauge Φ = exp(−R U(Ω)) Ψ, slider over R"),
    )
        defVar(ds, name, A, ("alpha", "theta", "R"); attrib = field_attrib(A,
            ["long_name" => ln, heat()..., slider...], angular))
    end
    logPsi = log10.(abs.(F["Psi"]))
    defVar(ds, "log10_abs_Psi", logPsi, ("alpha", "theta", "R"); attrib = field_attrib(logPsi,
        ["long_name" => "log10 |Ψ| (one colour scale over all R: exterior decay)", heat("Viridis")...,
         slider..., view_hint(playback_range = "sequence")...], angular))
    defVar(ds, "side", Int32.(F["side"]), ("R",); attrib = [
        "long_name" => "source of each R slice: −1 interior Fock series (R ≤ R*), +1 exterior Ritz (R > R*)",
        "flag_values" => Int32[-1, 1], "flag_meanings" => "interior exterior"])
    defVar(ds, "Rstar", fill(Rs), (); attrib = [
        "long_name" => "matching radius R*", "units" => "bohr", "H5WEB_geometry" => "vline",
        "H5WEB_label" => "R*", "H5WEB_color" => "#8a8985", "H5WEB_dash" => Int8(1)])
    # Ψ along R at three angle pairs (nearest grid points).
    picks = [(π / 2, π / 2), (π / 4, π / 2), (π / 2, 2π / 3)]
    idx = [(argmin(abs.(αp .- a)), argmin(abs.(θp .- t))) for (a, t) in picks]
    defVar(ds, "angle", collect(1:3), ("angle",); attrib = ["long_name" => "angle pair",
        "H5WEB_slice_labels" => "/angle_label"])
    defVar(ds, "angle_label", [@sprintf("(α, θ) = (%.3f, %.3f)", αp[i], θp[j]) for (i, j) in idx], ("angle",))
    along = hcat([F["Psi"][i, j, :] for (i, j) in idx]...)
    defVar(ds, "Psi_along_R", along, ("R", "angle"); attrib = [
        "long_name" => "Ψ along R at fixed angles (R* marked)",
        view_hint(vis = "line", scale = "log", x_dim = "R", multi_dim = "angle", overlays = "/Rstar")...])
    # Jumps across R* on the field grid.
    defVar(ds, "alpha_f", F["alpha"], ("alpha_f",); attrib = ["units" => "rad", "long_name" => "α", PI_TICKS])
    defVar(ds, "theta_f", F["theta"], ("theta_f",); attrib = ["units" => "rad", "long_name" => "θ", PI_TICKS])
    angular_f = write_coalescence!(ds; columns = "alpha_f theta_f", suffix = "_f")
    for (name, A, ln) in (
        ("jump_Psi_at_Rstar", M["Psi"] .- M["Psi_ext"], "Ψ(R*⁻) − Ψ(R*⁺): zero by construction (shared trace)"),
        ("jump_dRPsi_at_Rstar", M["dPsi_int"] .- M["dPsi_ext"], "∂_RΨ(R*⁻) − ∂_RΨ(R*⁺): the flux jump (matching residual)"),
    )
        defVar(ds, name, A, ("alpha_f", "theta_f"); attrib = field_attrib(A,
            ["long_name" => ln, sym...], angular_f; kind = :error))
    end
    # Interior series continued past R* minus the exterior solution.
    defVar(ds, "R_cont", F["R_cont"], ("R_cont",); attrib = ["units" => "bohr",
        "long_name" => "hyperradius R ≥ R* (interior series continued)"])
    C = F["cont_diff"]
    defVar(ds, "Psi_series_continued_minus_exterior", C, ("alpha", "theta", "R_cont"); attrib = field_attrib(C,
        ["long_name" => "Ψ_series(R) − Ψ_Ritz(R) for R ≥ R*: how far the two pieces agree off the sphere",
         sym..., view_hint(probe_dim = "R_cont", playback_dim = "R_cont")...], angular; kind = :error))
    maxC = [maximum(abs, C[:, :, k]) for k in axes(C, 3)]
    defVar(ds, "max_abs_series_minus_exterior", maxC, ("R_cont",); attrib = [
        "long_name" => "max over angles |Ψ_series − Ψ_Ritz| vs R ≥ R*", view_hint(vis = "line", scale = "log")...])
end
println("wrote ", path)
