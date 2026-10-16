# Recurrence residuals r_kp = (Λ² − k(k+4))ψ_kp − h_kp, with Λ² applied by ForwardDiff
# (not reconstructed from the recurrence), shown relative to the size of the terms,
# S = |Λ²ψ| + k(k+4)|ψ| + Σ|terms of h|, maximised over the grid (pointwise S vanishes on
# nodal lines of pure harmonics such as ψ21 ∝ cos θ):
#   * package evaluators, k ≤ 3 and ψ42 (ψ41 and ψ40 are Float64-only and not AD-able);
#   * spectral chain interpolant, every ψ_kp with k ≤ 8;
#   * package residual versus working precision (BigFloat) at three points.
# Writes output/hierarchy_residual.nc.
# Run: julia -t 4 --project=examples/sciviz examples/sciviz/hierarchy_residual.jl
include(joinpath(@__DIR__, "hierarchy_common.jl"))

const NG = 61
const KCHAIN = 8

"""Residual and scale at (α, θ) for coefficient (k, p), given value functions f(k, p)(α, θ)."""
function residual(fs, k, p, α, θ; Z, E)
    f = fs[(k, p)]
    val(kk, pp) = haskey(fs, (kk, pp)) ? fs[(kk, pp)](α, θ) : zero(α)
    L = Λ2(f, α, θ)
    v = f(α, θ)
    V = Vpot(α, θ; Z)
    h = recurrence_rhs(val, k, p, V, E)
    terms = (
        2(k + 2) * (p + 1) * val(k, p + 1),
        (p + 1) * (p + 2) * val(k, p + 2),
        2V * val(k - 1, p),
        2E * val(k - 2, p),
    )
    S = abs(L) + k * (k + 4) * abs(v) + sum(abs, terms)
    L - k * (k + 4) * v - h, S
end

function package_functions(T; Z, E, a21, rtol)
    Dict(
        (0, 0) => (a, t) -> one(a),
        (1, 0) => (a, t) -> psi10(a, t; Z),
        (2, 0) => (a, t) -> psi20(a, t; Z, E, a21),
        (2, 1) => (a, t) -> psi21(a, t; Z),
        (3, 0) => (a, t) -> psi30(a, t; Z, E, a21, rtol),
        (3, 1) => (a, t) -> psi31(a, t; Z),
        (4, 2) => (a, t) -> psi42(a, t; Z),
    )
end

αs = interior(NG; ε = 0.02)
θs = interior(NG; ε = 0.02)

# --- package, Float64 ---
pkg = package_functions(Float64; Z = HE.Z, E = HE.E, a21 = HE.a21, rtol = 1e-13)
pkeys = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (4, 2)]
pres = Dict{Tuple{Int,Int},Matrix{Float64}}()
for key in pkeys
    t = @elapsed begin
        RS = grid((a, b) -> residual(pkg, key..., a, b; Z = HE.Z, E = HE.E)[1], αs, θs)
        SS = grid((a, b) -> residual(pkg, key..., a, b; Z = HE.Z, E = HE.E)[2], αs, θs)
        pres[key] = abs.(RS) ./ maximum(SS)
    end
    @printf(
        "package ψ%d%d: median rel. residual %.1e, max %.1e (%.1f s)\n",
        key...,
        sort(vec(pres[key]))[end÷2],
        maximum(pres[key]),
        t
    )
end

# --- chain, Float64 ---
g, ψ = he_chain(; n = 64, kmax = KCHAIN)
x, bw = cheb_nodes(Float64, g.n)
ckeys = keys_upto(KCHAIN)
cfs = Dict(key => ((a, t) -> interp(x, bw, ψ[key], a, t)) for key in ckeys)
cres = Dict(key => Matrix{Float64}(undef, NG, NG) for key in ckeys if key != (0, 0))
cscale = Dict(key => Matrix{Float64}(undef, NG, NG) for key in ckeys if key != (0, 0))
t = @elapsed Threads.@threads for j = 1:NG
    for i = 1:NG, key in ckeys
        key == (0, 0) && continue
        r, S = residual(cfs, key..., αs[i], θs[j]; Z = HE.Z, E = HE.E)
        cres[key][i, j] = abs(r)
        cscale[key][i, j] = S
    end
end
for key in keys(cres)
    cres[key] ./= maximum(cscale[key])
end
@printf("chain residuals: %.1f s\n", t)
for key in ckeys
    key == (0, 0) && continue
    @printf(
        "chain ψ%d%d: median %.1e max %.1e\n",
        key...,
        sort(vec(cres[key]))[end÷2],
        maximum(cres[key])
    )
end

# --- package, BigFloat precision study ---
const PTS = [("1.1", "1.2"), ("0.4", "0.7"), ("2.2", "2.4")]
const BITS = [(53, "1e-13"), (128, "1e-30"), (192, "1e-48"), (256, "1e-66")]
bigres = Dict{Tuple{Int,Int},Vector{Float64}}(key => Float64[] for key in pkeys)
for (bits, tol) in BITS
    setprecision(BigFloat, max(bits, 64)) do
        T = bits == 53 ? Float64 : BigFloat
        fs = package_functions(
            T;
            Z = T(2),
            E = parse(T, "-2.903724377034119598311159245194404446696925310"),
            a21 = parse(T, "0.47674787900"),
            rtol = parse(T, tol),
        )
        for key in pkeys
            worst = 0.0
            for (a, b) in PTS
                r, S = residual(
                    fs,
                    key...,
                    parse(T, a),
                    parse(T, b);
                    Z = T(2),
                    E = parse(T, "-2.903724377034119598311159245194404446696925310"),
                )
                worst = max(worst, Float64(abs(r) / S))
            end
            push!(bigres[key], worst)
            @printf(
                "%3d bits ψ%d%d: max rel. residual over 3 points %.2e\n",
                bits,
                key...,
                worst
            )
        end
    end
end

# --- figures ---
crange = (-17, -9)
cmap = :magma
function resmap!(ax, A)
    heatmap!(ax, αs, θs, log10.(max.(A, 1e-18)); colormap = cmap, colorrange = crange)
end
ticks = ([0, π / 2, π], [L"0", L"\pi/2", L"\pi"])

fig = Figure(size = (1350, 760), fontsize = 15)
for (n, key) in enumerate(pkeys)
    r, c = fldmod1(n, 3)
    ax = Axis(
        fig[r, c];
        aspect = DataAspect(),
        xticks = ticks,
        yticks = ticks,
        title = latexstring("\\psi_{$(key[1]),$(key[2])}"),
        xlabel = r == 2 ? L"\alpha" : "",
        ylabel = c == 1 ? L"\theta" : "",
        xticklabelsvisible = r == 2,
        yticklabelsvisible = c == 1,
    )
    resmap!(ax, pres[key])
    if key == (2, 0)
        text!(
            ax,
            π / 2 + 0.08,
            2.6;
            text = "α = π/2:\nAD hits the\nremovable-limit\nbranch of χ",
            fontsize = 11,
            color = :white,
            align = (:left, :top),
        )
    end
end
Colorbar(
    fig[1:2, 4];
    colormap = cmap,
    colorrange = crange,
    label = L"\log_{10}\left(|r_{k,p}|/\max S_{k,p}\right)",
)
ax = Axis(
    fig[1:2, 5];
    yscale = log10,
    xlabel = L"\mathrm{working\ precision\ (bits)}",
    ylabel = L"\max\ |r_{k,p}|/S_{k,p}\ \ \mathrm{(3\ points,\ pointwise\ }S)",
    xticks = first.(BITS),
    title = L"\mathrm{BigFloat:\ residual\ vs\ precision}",
)
for (i, key) in enumerate(pkeys)
    nz = bigres[key] .> 0   # an exactly zero residual cannot be drawn on a log axis
    scatterlines!(
        ax,
        first.(BITS)[nz],
        bigres[key][nz];
        color = SERIES[i],
        linewidth = 2,
        markersize = 9,
        label = latexstring("\\psi_{$(key[1]),$(key[2])}"),
    )
end
lines!(
    ax,
    first.(BITS),
    2.0 .^ (-first.(BITS));
    color = :black,
    linestyle = :dot,
    label = L"2^{-\mathrm{bits}}",
)
axislegend(ax; position = :lb, framevisible = false, labelsize = 13)
colsize!(fig.layout, 5, Relative(0.3))
Label(
    fig[0, 1:3],
    L"\mathrm{FockExpansion.jl\ evaluators,\ ForwardDiff\ }\Lambda^2\mathrm{,\ Float64}";
    fontsize = 16,
)
savefig("hierarchy_residual_package", fig)

fig = Figure(size = (5 * 190 + 160, KCHAIN * 170 + 80), fontsize = 13)
for (k, p) in ckeys
    k == 0 && continue
    local ax = Axis(
        fig[k, p+1];
        aspect = DataAspect(),
        xticks = ticks,
        yticks = ticks,
        xticklabelsvisible = k == KCHAIN,
        yticklabelsvisible = p == 0,
        xticklabelsize = 11,
        yticklabelsize = 11,
        title = latexstring("\\psi_{$k,$p}"),
        titlegap = 2,
        titlesize = 13,
    )
    resmap!(ax, cres[(k, p)])
end
for k = 1:KCHAIN
    Label(fig[k, 0], latexstring("k=$k"); fontsize = 14, tellheight = false)
end
Colorbar(
    fig[1:3, 6];
    colormap = cmap,
    colorrange = crange,
    label = L"\log_{10}\left(|r_{k,p}|/\max S_{k,p}\right)",
    tellheight = false,
    height = Relative(0.9),
)
Label(
    fig[KCHAIN+1, 1:6],
    L"\mathrm{spectral\ chain\ }(n=64)\mathrm{\ interpolant;\ }\alpha\mathrm{\ horizontal,\ }\theta\mathrm{\ vertical}";
    fontsize = 13,
)
rowgap!(fig.layout, 4)
colgap!(fig.layout, 6)
savefig("hierarchy_residual_chain", fig)

# --- NetCDF ---
path = joinpath(OUT, "hierarchy_residual.nc")
nk, np = KCHAIN + 1, KCHAIN ÷ 2 + 1
C = fill(NaN, NG, NG, np, nk)
for (key, A) in cres
    C[:, :, key[2]+1, key[1]+1] = A
end
P = fill(NaN, NG, NG, length(pkeys))
for (i, key) in enumerate(pkeys)
    P[:, :, i] = pres[key]
end
isfile(path) && rm(path)
NCDataset(
    path,
    "c";
    attrib = [
        "title" => "Relative recurrence residuals |(Λ²−k(k+4))ψ_kp − h_kp| / S_kp (ForwardDiff)",
    ],
) do ds
    defVar(
        ds,
        "alpha",
        αs,
        ("alpha",);
        attrib = ["units" => "rad", "long_name" => "α (hyperangle)", PI_TICKS],
    )
    defVar(
        ds,
        "theta",
        θs,
        ("theta",);
        attrib = ["units" => "rad", "long_name" => "θ (interelectron angle)", PI_TICKS],
    )
    defVar(ds, "k", collect(0:KCHAIN), ("k",); attrib = ["long_name" => "power of R"])
    defVar(ds, "p", collect(0:(np-1)), ("p",); attrib = ["long_name" => "power of log R"])
    defVar(
        ds,
        "package_index",
        collect(1:length(pkeys)),
        ("package_index",);
        attrib = [
            "long_name" => "package evaluator",
            "comment" =>
                join(["$i: psi$(k)$(p)" for (i, (k, p)) in enumerate(pkeys)], ", "),
            "H5WEB_slice_labels" => "/package_name",
        ],
    )
    defVar(
        ds,
        "package_name",
        ["psi$(k)$(p)" for (k, p) in pkeys],
        ("package_index",);
        attrib = ["long_name" => "package evaluator name"],
    )
    hint = view_hint(vis = "heatmap", scale = "log", colormap = "Magma")
    angular = write_coalescence!(ds)
    defVar(
        ds,
        "chain_residual",
        C,
        ("alpha", "theta", "p", "k");
        attrib = field_attrib(
            C,
            [
                "long_name" => "relative residual, spectral chain n=64",
                hint...,
                view_hint(probe_dim = "k", playback_dim = "k")...,
            ],
            angular;
            kind = :error,
        ),
    )
    defVar(
        ds,
        "package_residual",
        P,
        ("alpha", "theta", "package_index");
        attrib = field_attrib(
            P,
            [
                "long_name" => "relative residual, package evaluators",
                hint...,
                view_hint(probe_dim = "package_index", playback_dim = "package_index")...,
            ],
            angular;
            kind = :error,
        ),
    )
    defVar(
        ds,
        "bits",
        first.(BITS),
        ("bits",);
        attrib = ["long_name" => "BigFloat precision (bits)"],
    )
    B = reduce(hcat, [bigres[key] for key in pkeys])
    defVar(
        ds,
        "bigfloat_residual",
        B,
        ("bits", "package_index");
        attrib = field_attrib(
            B,
            [
                "long_name" => "max relative residual over 3 points",
                view_hint(
                    vis = "line",
                    scale = "log",
                    envelope = "cummax",
                    envelope_direction = "reverse",
                )...,
            ];
            kind = :error,
        ),
    )
end
println("wrote ", path)
