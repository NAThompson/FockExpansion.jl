# steklov_energy.nc and steklov_constants.nc from data/energy.jls and data/constants.jls.
#     julia -t auto --project=examples/sciviz <this file>
include("/Users/nathompson7/FockExpansion.jl/examples/sciviz/common.jl")
using Printf, Serialization

const DATA = joinpath(@__DIR__, "data")
const STK = normpath(joinpath(@__DIR__, ".."))

# ---------------------------------------------------------------------------------------------
# steklov_energy.nc
# ---------------------------------------------------------------------------------------------
En = deserialize(joinpath(DATA, "energy.jls"))
dE, lam1, lam2, a21 = En["dE"], En["lambda_min"], En["lambda_2"], En["a21"]
# Least-squares line λ_min(E) = c0 + c1 dE and its zero.
X = [ones(length(dE)) dE]
c0, c1 = X \ lam1
dEstar = -c0 / c1
lam_fit = c0 .+ c1 .* dE
@printf("energy: dλ/dE = %.6e, zero crossing at E0 %+.3e, max fit residual %.2e\n", c1, dEstar,
        maximum(abs, lam1 .- lam_fit))
ca = (X \ a21)[2]

# The L = 10 three-point check of the study (energy_check.log).
chk = [(parse(Float64, m[1]), parse(Float64, m[2]), parse(Float64, m[3]), parse(Float64, m[4]))
       for m in (match(r"E = E0 ([-+\d.e]+): lambda_min = ([-+\d.e]+)\s+lambda_2 = ([-+\d.e]+)\s+a21 = ([\d.]+)", l)
                 for l in eachline(joinpath(STK, "energy_check.log"))) if m !== nothing]

path = joinpath(OUT, "steklov_energy.nc")
isfile(path) && rm(path)
setting = "R* = $(En["Rs"]), L = $(En["L"]), d = $(En["d"]), exterior $(En["basis"])"
NCDataset(path, "c"; attrib = [
    "title" => "Matching eigenvalues vs energy near E0",
    "setting" => setting, "E0" => En["E0"],
    "dlambda_dE" => c1, "zero_crossing_E_minus_E0" => dEstar, "da21_dE" => ca,
    "comment" => "λ_min > 0 below the crossing and < 0 above it. The crossing E* is a by-product energy (a diagnostic, not a bound).",
]) do ds
    defVar(ds, "dE", dE, ("dE",); attrib = ["units" => "hartree", "long_name" => "E − E0"])
    vl(name, x, label, color) = defVar(ds, name, fill(x), (); attrib = [
        "long_name" => label, "units" => "hartree", "H5WEB_geometry" => "vline",
        "H5WEB_label" => label, "H5WEB_color" => color, "H5WEB_dash" => Int8(1)])
    vl("E0_marker", 0.0, "E0", "#8a8985")
    vl("Estar_marker", dEstar, @sprintf("λ_min = 0 at E0 %+.2e", dEstar), "#eb6834")
    defVar(ds, "lambda_zero", fill(0.0), (); attrib = ["long_name" => "λ = 0",
        "H5WEB_geometry" => "hline", "H5WEB_label" => "λ = 0", "H5WEB_color" => "#8a8985"])
    defVar(ds, "lambda_min_fit", lam_fit, ("dE",); attrib = [
        "long_name" => @sprintf("linear fit, slope %.4e", c1), "H5WEB_label" => "linear fit",
        "H5WEB_color" => "#8a8985", "H5WEB_dash" => Int8(1)])
    defVar(ds, "lambda_min", lam1, ("dE",); attrib = [
        "long_name" => "smallest-|λ| eigenvalue of M⁻¹(B_int + A_ext) vs E ($setting)", "units" => "1",
        view_hint(vis = "line", series = "/lambda_min_fit", overlays = "/E0_marker /Estar_marker /lambda_zero")...])
    defVar(ds, "lambda_min_minus_fit", lam1 .- lam_fit, ("dE",); attrib = [
        "long_name" => "λ_min − linear fit (curvature / noise)", view_hint(vis = "line")...])
    defVar(ds, "lambda_2", lam2, ("dE",); attrib = [
        "long_name" => "second-smallest |λ| (gap: the matched mode is isolated)", "units" => "1",
        view_hint(vis = "line", overlays = "/E0_marker")...])
    defVar(ds, "a21", a21, ("dE",); attrib = [
        "long_name" => @sprintf("a21 from the λ_min eigenvector vs E (slope %.4f per hartree)", ca),
        view_hint(vis = "line", overlays = "/E0_marker")...])
    defVar(ds, "seconds", En["seconds"], ("dE",); attrib = ["long_name" => "wall time per point", "units" => "s"])
    g = defGroup(ds, "L10_check")
    defVar(g, "dE", [c[1] for c in chk], ("dE",); attrib = ["units" => "hartree", "long_name" => "E − E0"])
    defVar(g, "lambda_min", [c[2] for c in chk], ("dE",); attrib = [
        "long_name" => "λ_min, study check L = 10, d = 32, B10P6 (energy_check.log)", view_hint(vis = "line")...])
    defVar(g, "lambda_2", [c[3] for c in chk], ("dE",); attrib = ["long_name" => "λ_2, L = 10 check"])
    defVar(g, "a21", [c[4] for c in chk], ("dE",); attrib = ["long_name" => "a21, L = 10 check", view_hint(vis = "line")...])
end
println("wrote ", path)

# ---------------------------------------------------------------------------------------------
# steklov_constants.nc
# ---------------------------------------------------------------------------------------------
C = deserialize(joinpath(DATA, "constants.jls"))
names = C["names"]
nc = length(names)
row(k) = C[k]
p_err = 4.87
tail(x44, x48) = (x48 .- x44) ./ ((48 / 44)^p_err - 1)   # study's algebraic d-tail model
v = Dict(k => row(k) for k in keys(C) if occursin("d=", k))
best48 = v["L=14 B10P6 d=48"]
best = best48 .+ tail(v["L=14 B10P6 d=44"], best48)
dtail = best .- best48
Lstep = v["L=14 B10P6 d=48"] .- v["L=12 B10P6 d=48"]
bstep = v["L=12 B12P6 d=48"] .- v["L=12 B10P6 d=48"]

# R* spread at L = 12, d = 40, B10P6 (logs): a21, a_(2,0), a_(0,2) only.
function log_row(file)
    for l in eachline(joinpath(STK, file))
        m = match(r"L=12 B10P6 d=40: .*a21=([\d.]+) a3=([-+\d.]+) a4=([-+\d.]+)", l)
        m === nothing || return parse.(Float64, (m[1], m[2], m[3]))
    end
end
rows = [log_row(f) for f in ("sweep_d_R075_L12.log", "sweep_d_R1_B10P6_L12.log", "sweep_d_R15_L12.log")]
rspread = fill(NaN, nc)
for (i, col) in ((1, 1), (3, 2), (4, 3))
    x = [r[col] for r in rows]
    rspread[i] = maximum(x) - minimum(x)
end
unc = [maximum(filter(isfinite, abs.([dtail[i], Lstep[i], bstep[i], rspread[i]]))) for i = 1:nc]
decimals = [floor(Int, -log10(u)) for u in unc]
show_value(x, u) = (n = floor(Int, -log10(u)); s = @sprintf("%.*f", n + 1, x); "$s ± $(@sprintf("%.0e", u))")
summary = ["$(names[i]) = $(show_value(best[i], unc[i]))" for i = 1:nc]
foreach(println, summary)

path = joinpath(OUT, "steklov_constants.nc")
isfile(path) && rm(path)
NCDataset(path, "c"; attrib = [
    "title" => "Fock constants of helium from Steklov–Poincaré matching (R* = 1)",
    "best_setting" => "L = 14, B10P6, d = 48 plus the algebraic d-tail (error ∝ d^-4.87)",
    "reference_a21" => 0.476747879,
    "comment" => "Cusp-gauge seeds a_(m,j): coefficient of R^K Y_(m,j) in Φ = exp(−S)Ψ, K = 2(m+j), a_(0,0) = 1. a21 and <psi40, Y> are physical total projections. uncertainty = max(|d tail|, |L 12→14|, |basis B12P6−B10P6|, R* spread where known); provisional, not a bound.",
]) do ds
    defVar(ds, "constant", collect(1:nc), ("constant",); attrib = [
        "long_name" => "constant index", "H5WEB_slice_labels" => "/label"])
    defVar(ds, "label", names, ("constant",); attrib = ["long_name" => "constant name"])
    defVar(ds, "order_K", [2, 2, 4, 4, 6, 6, 8, 8, 8, 4, 4], ("constant",); attrib = [
        "long_name" => "Fock order K (power of R)"])
    defVar(ds, "summary", summary, ("constant",); attrib = ["long_name" => "value ± uncertainty (stable digits shown + 1)"])
    put(name, x, ln) = defVar(ds, name, x, ("constant",); attrib = ["long_name" => ln])
    put("value_best", best, "best estimate: L = 14 B10P6 d = 48 + d-tail")
    put("uncertainty", unc, "provisional uncertainty (max of the components below)")
    defVar(ds, "stable_decimals", Int32.(decimals), ("constant",); attrib = [
        "long_name" => "number of stable decimal places, floor(−log10 uncertainty)"])
    put("d_tail", dtail, "algebraic d-tail added to d = 48 (L = 14 B10P6)")
    put("L_step_12_to_14", Lstep, "a(L=14) − a(L=12), B10P6, d = 48")
    put("basis_B12P6_minus_B10P6", bstep, "exterior basis step at L = 12, d = 48")
    put("Rstar_spread", rspread, "max − min over R* ∈ {0.75, 1, 1.5} at L = 12, d = 40 (NaN = not logged)")
    put("rel_uncertainty", unc ./ abs.(best), "uncertainty / |value|")
    for k in sort(collect(keys(v)))
        put("value_" * replace(k, " " => "_", "=" => ""), v[k], "value at $k (R* = 1)")
    end
end
println("wrote ", path)
