# Shared helpers for the paper figures. Run scripts with
#     julia -t auto --project=examples/sciviz examples/sciviz/<script>.jl
# Each script writes PDF/PNG figures for the paper and a NetCDF-4 (= HDF5) file for the
# H5Web viewer (see VISUALIZATION_NOTES.md) into examples/sciviz/output.
using FockExpansion, CairoMakie, NCDatasets

const OUT = joinpath(@__DIR__, "output")
mkpath(OUT)

# Helium ground-state inputs (README, "ForwardDiff recurrence residual").
const HE = (Z = 2.0, E = -2.903724377034119598311159245194404446696925310, a21 = 0.47674787900)

"""Interior grid on (0, π) with `n` points, avoiding the coordinate singularities."""
interior(n; ε = 1e-3) = collect(range(ε, π - ε; length = n))

"""Evaluate `f(α, θ)` on the grid, threaded over θ. Returns a matrix indexed [iα, iθ]."""
function grid(f, αs, θs)
    A = Matrix{Float64}(undef, length(αs), length(θs))
    Threads.@threads for j in eachindex(θs)
        for i in eachindex(αs)
            A[i, j] = f(αs[i], θs[j])
        end
    end
    A
end

"""Geometry of the electron triangle: r1/ρ = cos(α/2), r2/ρ = sin(α/2),
r12/ρ = sqrt(1 - sin α cos θ)."""
hypersphere_xyz(α, θ) = (cos(α / 2), sin(α / 2) * cos(θ), sin(α / 2) * sin(θ))

"""Typed H5Web view attributes (v1 schema, plain HDF5 attributes):
`view_hint(vis = "heatmap", scale = "log")` gives
`["H5WEB_vis" => "heatmap", "H5WEB_scale" => "log"]`; splat it into an `attrib` vector.
Keywords are the attribute names without the prefix (`vis`, `scale`, `colormap`,
`domain`, `symmetric`, `center`, `keep_ratio`, `overlays`, ...). Booleans are written as
Int8 (NetCDF has no Bool), `layout` as Int32 and other vectors (`domain = [lo, hi]`) as
Float64, with NaN or `nothing` for an automatic end."""
view_hint(; kw...) = Pair{String,Any}["H5WEB_$k" => h5web_value(k, v) for (k, v) in kw]

h5web_value(k, v::Bool) = Int8(v)
h5web_value(k, v::AbstractVector) =
    k === :layout ? Int32.(v) : Float64[something(x, NaN) for x in v]
h5web_value(k, v) = v

# Coalescence loci, drawn as overlays on every angular field: the electron-nucleus lines
# r1 = 0 (α = π) and r2 = 0 (α = 0) and the electron-electron point (π/2, 0).
const COALESCENCE = [
    ("electron_nucleus_1", "curve", "r_{1} = 0", [π 0.0; π π]),
    ("electron_nucleus_2", "curve", "r_{2} = 0", [0.0 0.0; 0.0 π]),
    ("electron_electron", "points", "r_{12} = 0", [π/2 0.0]),
]

# The same loci in the spectral chain's separating coordinates X = α − β, Y = α + β
# (cos β = sin α cos θ), where all three are corners of the grid.
const COALESCENCE_XY = [
    ("electron_nucleus_1_XY", "points", "r_{1} = 0", [π/2 3π/2]),
    ("electron_nucleus_2_XY", "points", "r_{2} = 0", [-π/2 π/2]),
    ("electron_electron_XY", "points", "r_{12} = 0", [π/2 π/2]),
]

"""Axis attribute: ticks at multiples of π on an angle dimension scale."""
const PI_TICKS = "H5WEB_ticks" => "pi"

"""Write `loci` (as `COALESCENCE`, names + `suffix`) into the group `coalescence` of `ds`
as (N, 2) overlay datasets with their own style, `columns` naming the base's dimension
scales; existing ones are kept. `prefix` is the path of `ds` when it is a subgroup
(e.g. "/growth"). Returns the default view attributes of a field they are drawn on
(automatic colour map, equal aspect, these overlays); combine them with the field's own
via `with_defaults`."""
function write_coalescence!(
    ds,
    loci = COALESCENCE;
    columns = "alpha theta",
    suffix = "",
    prefix = "",
)
    coal =
        haskey(ds.group, "coalescence") ? ds.group["coalescence"] :
        defGroup(ds, "coalescence")
    haskey(coal.dim, "pair") || defDim(coal, "pair", 2)
    names = [c[1] * suffix for c in loci]
    for (cname, (_, geometry, label, pts)) in zip(names, loci)
        haskey(coal, cname) && continue
        defDim(coal, cname * "_n", size(pts, 1))
        # Julia is column-major: a (2, N) array appears in the file as (N, 2).
        defVar(
            coal,
            cname,
            permutedims(pts),
            ("pair", cname * "_n");
            attrib = [
                "long_name" => "coalescence $label",
                "H5WEB_geometry" => geometry,
                "H5WEB_columns" => columns,
                "H5WEB_color" => "#ff2fd6",
                # Lines take a width, points a radius.
                geometry == "points" ? ("H5WEB_size" => 8.0) : ("H5WEB_width" => 2.0),
                "H5WEB_label" => label,
                "H5WEB_label_mode" => "inline",
            ],
        )
    end
    view_hint(
        colormap = "auto",
        keep_ratio = true,
        overlays = join(["$prefix/coalescence/$n" for n in names], " "),
    )
end

"""Names of residuals, errors and differences (`kind = :error` by default)."""
const ERROR_NAME = r"residual|err|diff|minus|delta"

field_kind(name) = occursin(ERROR_NAME, name) ? :error : :field

"""View attributes read off the data `A` (NaN ignored), for a field whose own attributes
are `own`. `kind = :field` (smooth physical fields): the nodal line emphasized
(`H5WEB_contour_emphasize = [0]`) when `A` changes sign, and a symlog threshold at the
power of 10 nearest 1e-3 times the 99th percentile of |A| (robust to coalescence
spikes). `kind = :error` (residuals, errors, differences): no emphasized contour (it would trace noise), and a symlog threshold at the power of 10
nearest the 1st percentile of the nonzero |A|, clamped to [1e-16, 1e-2]. Either kind:
a grey `H5WEB_nan_color` when `A` has NaN (or, on a log scale, zeros), and `H5WEB_playback_skip_empty` when one of
its images (the first two dimensions, at fixed indices of the rest) is entirely NaN."""
function data_hints(A, own = Pair{String,Any}[]; kind = :field)
    kind in (:field, :error) || throw(ArgumentError("kind = $kind"))
    v = filter(isfinite, vec(A))
    hints = Pair{String,Any}[]
    scale = [last(a) for a in own if first(a) == "H5WEB_scale"]
    # Missing values, and zeros on a log scale, are drawn grey.
    if length(v) < length(A) || ("log" in scale && any(iszero, v))
        append!(hints, view_hint(nan_color = "#d9d9d9"))
    end
    if length(v) < length(A)
        images = reshape(A, size(A, 1) * size(A, 2), :)
        ndims(A) > 2 && any(j -> all(isnan, view(images, :, j)), axes(images, 2)) &&
            append!(hints, view_hint(playback_skip_empty = true))
    end
    isempty(v) && return hints
    lo, hi = extrema(v)
    kind == :field && lo < 0 < hi &&
        append!(hints, view_hint(contour_emphasize = [0.0]))
    if "symlog" in scale
        m = sort!(filter!(!iszero, abs.(v)))
        if isempty(m)
            t = 1.0
        elseif kind == :field
            t = exp10(round(log10(1e-3 * m[ceil(Int, 0.99 * length(m))])))
        else
            q = m[max(1, ceil(Int, 0.01 * length(m)))]
            t = clamp(exp10(round(log10(q))), 1e-16, 1e-2)
        end
        append!(hints, view_hint(symlog_threshold = t))
    end
    hints
end

"""Axis attributes of a dimension scale `x`: `H5WEB_axis_scale = "log"` when it is
positive and geometric (ratios of neighbours equal to within 1%)."""
function axis_hint(x)
    r = x[2:end] ./ x[1:end-1]
    geometric = length(x) >= 3 && all(>(0), x) && maximum(r) <= 1.01 * minimum(r)
    geometric ? view_hint(axis_scale = "log") : Pair{String,Any}[]
end

"""`own` attributes followed by the `defaults` whose names `own` does not set."""
with_defaults(own, defaults) =
    Pair{String,Any}[own; filter(d -> first(d) ∉ first.(own), defaults)]

"""Attributes of a field with data `A`: its `own`, then `data_hints(A, own; kind)` and
the `defaults` that `own` does not override."""
field_attrib(A, own, defaults = Pair{String,Any}[]; kind = :field) =
    with_defaults(own, [data_hints(A, own; kind); defaults])

"""Write the hypersphere coordinates of the (αs, θs) grid with dimensions `dims` into `ds`
as `x`, `y`, `z` + `suffix`; returns their names for a CF `coordinates` attribute."""
function write_xyz!(ds, αs, θs, dims; suffix = "")
    long_names = ("x = r_{1}/ρ", "y = (r_{2}/ρ) cos θ", "z = (r_{2}/ρ) sin θ")
    names = ["x", "y", "z"] .* suffix
    for (c, (name, long_name)) in enumerate(zip(names, long_names))
        X = [hypersphere_xyz(a, t)[c] for a in αs, t in θs]
        defVar(ds, name, X, dims; attrib = ["long_name" => long_name])
    end
    join(names, " ")
end

"""Write angular fields `fields::Vector{Pair{String,Matrix}}` (each [iα, iθ]) to a NetCDF
file with `alpha`, `theta` dimension scales and the hypersphere coordinates `x y z`, so
each field opens as a heatmap and also as a parametric surface. `attribs` maps a field
name to extra attributes, e.g. `view_hint(vis = "heatmap")`; they override the default
view attributes (automatic colour map, equal aspect, coalescence overlays)."""
function write_angular(
    path,
    αs,
    θs,
    fields;
    attribs = Dict{String,Any}(),
    global_attribs = Pair{String,Any}[],
)
    isfile(path) && rm(path)
    NCDataset(path, "c"; attrib = global_attribs) do ds
        defVar(
            ds,
            "alpha",
            αs,
            ("alpha",);
            attrib = [
                "units" => "rad",
                "long_name" => "α (hyperangle)",
                PI_TICKS,
            ],
        )
        defVar(
            ds,
            "theta",
            θs,
            ("theta",);
            attrib = [
                "units" => "rad",
                "long_name" => "θ (angle between r1, r2)",
                PI_TICKS,
            ],
        )
        write_xyz!(ds, αs, θs, ("alpha", "theta"))
        defaults = write_coalescence!(ds)
        for (name, A) in fields
            own = get(attribs, name, Pair{String,Any}[])
            own = ["coordinates" => "x y z"; own]
            at = field_attrib(A, own, defaults; kind = field_kind(name))
            defVar(ds, name, A, ("alpha", "theta"); attrib = at)
        end
    end
    println("wrote ", path)
end

function savefig(name, fig)
    for ext in ("pdf", "png")
        save(joinpath(OUT, name * "." * ext), fig; px_per_unit = 2)
    end
    println("saved ", name)
end

set_theme!(theme_latexfonts(); fontsize = 16)
