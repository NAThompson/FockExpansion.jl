# steklov_convergence.nc: a21 convergence of the Steklov–Poincaré matching, from the study logs.
#     julia -t auto --project=examples/sciviz <this file>
include("/Users/nathompson7/FockExpansion.jl/examples/sciviz/common.jl")
using Printf

const STK = normpath(joinpath(@__DIR__, ".."))
const REF = 0.476747879
const STEKLOV_VALUE = 0.47674788   # study's reported value, provisional ±1e-8

# Codex fixed-energy hyperspherical-channel sequence (r = 0.75), d = 4..24.
const CODEX_D = [4, 8, 12, 16, 20, 24]
const CODEX_A21 = [0.476594728101, 0.476737464886, 0.476746090427, 0.476747395497,
                   0.476747707937, 0.476747806703]

# "R*=1.00 L=12 B12P6 [nb= 894] d=24: lambda=+1.6e-08 a21=0.47..." (q=… quadrature lines skipped).
function parse_log(file)
    rows = NamedTuple[]
    for line in eachline(joinpath(STK, file))
        m = match(
            r"R\*=([\d.]+) L=\s*(\d+) (\S+)\s+(?:nb=\s*\d+ )?d=\s*(\d+): lambda=([-+\d.e]+) a21=([\d.]+)",
            line,
        )
        m === nothing && continue
        push!(rows, (Rs = parse(Float64, m[1]), L = parse(Int, m[2]), basis = m[3],
                     d = parse(Int, m[4]), lambda = parse(Float64, m[5]), a21 = parse(Float64, m[6])))
    end
    return rows
end

const LOGS = ["sweep_d_R1.log", "sweep_d_R075.log", "sweep_d_R15.log", "sweep_d_R05.log",
              "sweep_d_R2.log", "sweep_d_R075_L12.log", "sweep_d_R15_L12.log",
              "sweep_d_R1_B10P6.log", "sweep_d_R1_B10P6_L12.log", "extbasis_R1.log",
              "extbasis_R1_hiL.log", "extbasis_R1_L14.log"]
rows = unique(r -> (r.Rs, r.L, r.basis, r.d), vcat(parse_log.(LOGS)...))
println(length(rows), " distinct (R*, L, basis, d) runs parsed")

lookup(Rs, L, basis, d) = (i = findfirst(r -> (r.Rs, r.L, r.basis, r.d) == (Rs, L, basis, d), rows);
                           i === nothing ? nothing : rows[i])
value(Rs, L, basis, d, f = :a21) = (r = lookup(Rs, L, basis, d); r === nothing ? NaN : getfield(r, f))
fmtR(Rs) = rstrip(rstrip(@sprintf("%.2f", Rs), '0'), '.')

# ---------------------------------------------------------------------------------------------
# Curve families
# ---------------------------------------------------------------------------------------------
dgrid = collect(4:4:48)
# (label, R*, L, basis); the Codex row is filled separately.
d_series = [(Rs, L, b) for (Rs, L, b) in sort(unique((r.Rs, r.L, r.basis) for r in rows))
            if count(r -> (r.Rs, r.L, r.basis) == (Rs, L, b), rows) >= 3]
d_labels = vcat("Codex hyperspherical channels (r = 0.75)",
                ["Steklov R* = $(fmtR(Rs)), L = $L, $b" for (Rs, L, b) in d_series])
A_d = fill(NaN, length(dgrid), length(d_labels))      # Julia (d, series) → HDF5 (series, d)
Lam_d = fill(NaN, size(A_d))
for (i, d) in enumerate(dgrid)
    j = findfirst(==(d), CODEX_D)
    j === nothing || (A_d[i, 1] = CODEX_A21[j])
    for (k, (Rs, L, b)) in enumerate(d_series)
        A_d[i, k+1] = value(Rs, L, b, d)
        Lam_d[i, k+1] = value(Rs, L, b, d, :lambda)
    end
end

# Increments per unit d, |a21(d) − a21(d_prev)| / (d − d_prev), placed at d.
Inc = fill(NaN, size(A_d))
for k in axes(A_d, 2)
    idx = findall(isfinite, A_d[:, k])
    for (p, q) in zip(idx[1:end-1], idx[2:end])
        Inc[q, k] = abs(A_d[q, k] - A_d[p, k]) / (dgrid[q] - dgrid[p])
    end
end
# Guide ∝ d^-6.87 per unit d (= the study's d^-5.87 law for steps of 4, divided by 4).
guide = [d >= 12 ? 2.5e-7 / 4 * (d / 18)^-5.87 : NaN for d in dgrid]

# a21 vs trace degree L.
Lgrid = collect(4:2:14)
L_series = [(1.0, "B12P6", 24), (1.0, "B12P6", 48), (1.0, "B10P6", 24), (1.0, "B10P6", 40),
            (1.0, "B10P6", 48), (0.75, "B10P6", 24), (0.75, "B10P6", 40), (1.5, "B10P6", 24),
            (1.5, "B10P6", 40), (0.5, "B10P6", 32), (2.0, "B10P6", 32)]
L_labels = ["R* = $(fmtR(Rs)), d = $d, $b" for (Rs, b, d) in L_series]
A_L = [value(Rs, L, b, d) for L in Lgrid, (Rs, b, d) in L_series]
Lam_L = [value(Rs, L, b, d, :lambda) for L in Lgrid, (Rs, b, d) in L_series]

# Exterior basis spread at R* = 1: a21(basis) − a21(B12P6) for d = 16 and d = 24.
bases = ["B8P2", "B10P4", "B10P6", "B12P4", "B12P6"]
spread_labels = ["$b, d = $d" for d in (16, 24) for b in bases]
spread_keys = [(b, d) for d in (16, 24) for b in bases]
Spread = [value(1.0, L, b, d) - value(1.0, L, "B12P6", d) for L in Lgrid, (b, d) in spread_keys]
LamSpread = [value(1.0, L, b, d, :lambda) for L in Lgrid, (b, d) in spread_keys]

# a21 vs matching radius R*.
Rgrid = [0.5, 0.75, 1.0, 1.5, 2.0]
R_series = [(8, "B10P6", 32), (10, "B10P6", 32), (10, "B10P6", 40), (12, "B10P6", 40)]
R_labels = ["L = $L, d = $d, $b" for (L, b, d) in R_series]
A_R = [value(Rs, L, b, d) for Rs in Rgrid, (L, b, d) in R_series]

# ---------------------------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------------------------
path = joinpath(OUT, "steklov_convergence.nc")
isfile(path) && rm(path)

# A family dimension with string labels used as slice labels / legend.
function family!(g, prefix, dim, labels, comment)
    defVar(g, dim, collect(1:length(labels)), (dim,);
           attrib = ["long_name" => comment, "H5WEB_slice_labels" => "$prefix/$(dim)_label"])
    defVar(g, "$(dim)_label", labels, (dim,); attrib = ["long_name" => "curve labels"])
end

# Horizontal reference line (scalar dataset).
function hline!(g, name, y, label; color = "#52514e", dash = true)
    defVar(g, name, fill(Float64(y)), (); attrib = [
        "long_name" => label, "H5WEB_geometry" => "hline", "H5WEB_label" => label,
        "H5WEB_color" => color, "H5WEB_dash" => Int8(dash)])
end

function yband!(g, name, lo, hi, label; color = "#eda100")
    defDim(g, "$(name)_end", 2)
    defVar(g, name, [lo, hi], ("$(name)_end",); attrib = [
        "long_name" => "$label [lo, hi]", "H5WEB_geometry" => "band", "H5WEB_axis" => "y",
        "H5WEB_label" => label, "H5WEB_color" => color])
end

lineview(prefix, x, multi; scale = "linear", overlays = String[]) = view_hint(
    vis = "line", scale = scale, x_dim = x, multi_dim = multi,
    overlays = join(["$prefix/$o" for o in overlays], " "))

NCDataset(path, "c"; attrib = [
    "title" => "Steklov–Poincaré matching for helium at E0: a21 convergence",
    "source" => "parsed from $(join(LOGS, ", ")) in $STK; Codex channel sequence from figure.jl",
    "reference_a21" => REF,
    "E0" => -2.903724377034119598311,
    "comment" => "a21 = <psi20, s>/<s, s> (physical total projection); lambda = smallest-|.| eigenvalue of M^-1 (B_int + A_ext)",
]) do ds
    # --- /by_d
    g = defGroup(ds, "by_d")
    defVar(g, "d", dgrid, ("d",); attrib = ["long_name" => "interior channel degree d", "units" => "1"])
    defVar(g, "channels", [(d + 2)^2 ÷ 4 for d in dgrid], ("d",);
           attrib = ["long_name" => "number of interior channels ≈ (d+2)²/4"])
    family!(g, "/by_d", "series", d_labels, "curve (Codex channels or Steklov R*, L, exterior basis)")
    hline!(g, "reference", REF, "reference a21 = 0.476747879")
    yband!(g, "steklov_estimate", STEKLOV_VALUE - 1e-8, STEKLOV_VALUE + 1e-8, "study estimate 0.47674788 ± 1e-8")
    hline!(g, "zero", 0.0, "reference (difference = 0)")
    defVar(g, "a21", A_d, ("d", "series"); attrib = [
        "long_name" => "a21 vs interior degree d (NaN = not run)", "units" => "1",
        lineview("/by_d", "d", "series"; overlays = ["reference", "steklov_estimate"])...])
    defVar(g, "a21_minus_reference", A_d .- REF, ("d", "series"); attrib = [
        "long_name" => "a21 − 0.476747879 vs interior degree d", "units" => "1",
        lineview("/by_d", "d", "series"; scale = "symlog", overlays = ["zero"])...,
        "H5WEB_symlog_threshold" => 1e-9])
    defVar(g, "lambda_min", Lam_d, ("d", "series"); attrib = [
        "long_name" => "smallest matching eigenvalue λ_min (combined matching residual; Codex rows NaN)",
        "units" => "1",
        lineview("/by_d", "d", "series"; scale = "log")...])

    # --- /increments (log–log)
    g = defGroup(ds, "increments")
    defVar(g, "d", dgrid, ("d",); attrib = ["long_name" => "interior channel degree d",
                                            view_hint(axis_scale = "log")...])
    family!(g, "/increments", "series", d_labels, "curve")
    defVar(g, "guide_d_power", guide, ("d",); attrib = [
        "long_name" => "guide ∝ d^-6.87 per unit d (study's d^-5.87 per step of 4)",
        "H5WEB_label" => "∝ d^{-6.87}", "H5WEB_color" => "#8a8985", "H5WEB_dash" => Int8(1)])
    defVar(g, "abs_increment_per_unit_d", Inc, ("d", "series"); attrib = [
        "long_name" => "|a21(d) − a21(d_prev)| / (d − d_prev), at d (steps of 4 or 8)",
        "units" => "1",
        lineview("/increments", "d", "series"; scale = "log")...])
    defVar(g, "pair", [1, 2], ("pair",); attrib = ["long_name" => "1 Codex, 2 Steklov L=12",
                                                     "H5WEB_slice_labels" => "/increments/pair_label"])
    defVar(g, "pair_label", ["Codex channels", "Steklov R* = 1, L = 12, B12P6"], ("pair",))
    defVar(g, "abs_increment_main", Inc[:, [1; findfirst(==("Steklov R* = 1, L = 12, B12P6"), d_labels)]],
           ("d", "pair"); attrib = [
        "long_name" => "|Δa21|/Δd: Codex channels vs Steklov R* = 1, L = 12, B12P6 (with d^-6.87 guide)",
        view_hint(vis = "line", scale = "log", x_dim = "d", multi_dim = "pair",
                  series = "/increments/guide_d_power")...])

    # --- /by_L
    g = defGroup(ds, "by_L")
    defVar(g, "L", Lgrid, ("L",); attrib = ["long_name" => "trace harmonic degree L (K ≤ 2L)"])
    family!(g, "/by_L", "series", L_labels, "curve (R*, interior d, exterior basis)")
    hline!(g, "reference", REF, "reference a21 = 0.476747879")
    yband!(g, "steklov_estimate", STEKLOV_VALUE - 1e-8, STEKLOV_VALUE + 1e-8, "study estimate 0.47674788 ± 1e-8")
    defVar(g, "a21", A_L, ("L", "series"); attrib = [
        "long_name" => "a21 vs trace degree L (NaN = not run)", "units" => "1",
        lineview("/by_L", "L", "series"; overlays = ["reference", "steklov_estimate"])...])
    defVar(g, "lambda_min", Lam_L, ("L", "series"); attrib = [
        "long_name" => "λ_min vs trace degree L", "units" => "1",
        lineview("/by_L", "L", "series"; scale = "log")...])

    # --- /exterior_basis
    g = defGroup(ds, "exterior_basis")
    defVar(g, "L", Lgrid, ("L",); attrib = ["long_name" => "trace harmonic degree L"])
    family!(g, "/exterior_basis", "series", spread_labels, "exterior basis BxPy: Hylleraas Ω = x, y channel corrections")
    hline!(g, "zero", 0.0, "B12P6 (largest basis)")
    defVar(g, "a21_minus_B12P6", Spread, ("L", "series"); attrib = [
        "long_name" => "a21(basis) − a21(B12P6), R* = 1, same d and L", "units" => "1",
        lineview("/exterior_basis", "L", "series"; scale = "symlog", overlays = ["zero"])...,
        "H5WEB_symlog_threshold" => 1e-10])
    defVar(g, "lambda_min", LamSpread, ("L", "series"); attrib = [
        "long_name" => "λ_min per exterior basis, R* = 1", "units" => "1",
        lineview("/exterior_basis", "L", "series"; scale = "log")...])

    # --- /by_Rstar
    g = defGroup(ds, "by_Rstar")
    defVar(g, "Rstar", Rgrid, ("Rstar",); attrib = ["long_name" => "matching radius R*", "units" => "bohr"])
    family!(g, "/by_Rstar", "series", R_labels, "curve (L, d, exterior basis)")
    hline!(g, "reference", REF, "reference a21 = 0.476747879")
    defVar(g, "a21", A_R, ("Rstar", "series"); attrib = [
        "long_name" => "a21 vs matching radius R* (NaN = not run)", "units" => "1",
        lineview("/by_Rstar", "Rstar", "series"; overlays = ["reference"])...])
end
println("wrote ", path)
