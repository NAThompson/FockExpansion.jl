# Schematic of Codex's constructive patched representation (EFFECTIVE_* notes):
# nested radial regions on a log2 ρ axis, and the pipeline
#   certified trial → Temple → heat smoothing → clip/Dyson/cubature → W
#   certified trial → averaged trace inversion → recurrence → F_N → blend.
#
#     julia -t 4 --project=examples/sciviz examples/sciviz/patched_schematic.jl
include(joinpath(@__DIR__, "patched_common.jl"))

seed = read(joinpath(BM, "EFFECTIVE_SEED_CERTIFICATE.json"), String)
heat = read(joinpath(BM, "EFFECTIVE_HEAT_BUDGET.json"), String)
Ms = parse.(Int, jsonall(heat, "clipping_M"))
epsk = [round(Int, log2(parse(Int, split(s, "/")[2]))) for s in jsonall(heat, "epsilon")]

l2(x) = log2(Float64(x))
xr = (-22.5, 3.5)
fig = Figure(size = (1500, 1150))

# ------------------------------------------------------------------ weights
axw = Axis(
    fig[1, 1];
    height = 110,
    xticksvisible = false,
    xticklabelsvisible = false,
    ylabel = "weight",
    yticks = [0, 1],
    title = L"(a) radial regions of the patched approximant $Z_k = \chi F_{N,k} + (1-\chi) W_k$ ($\rho$ in bohr)",
    titlesize = 17,
)
xs = range(xr...; length = 2000)
lines!(
    axw,
    xs,
    [chi(2.0^x, 2.0^-14) for x in xs];
    color = C.blue,
    linewidth = 2.5,
    label = L"\chi",
)
lines!(
    axw,
    xs,
    [1 - chi(2.0^x, 2.0^-14) for x in xs];
    color = C.gray,
    linewidth = 2.5,
    label = L"1-\chi",
)
axislegend(axw; position = :rc, framevisible = false, orientation = :horizontal)
xlims!(axw, xr...)

# ------------------------------------------------------------------ regions
axr = Axis(
    fig[2, 1];
    height = 380,
    xlabel = L"\log_2 \rho",
    yticksvisible = false,
    yticklabelsvisible = false,
    xticks = -22:2:2,
    ygridvisible = false,
)
linkxaxes!(axw, axr)
xlims!(axr, xr...)
ylims!(axr, -0.6, 6.0)
function bar!(ax, x0, x1, y, col; alpha = 0.85, h = 0.32)
    poly!(ax, Rect(x0, y - h, x1 - x0, 2h); color = (col, alpha), strokewidth = 0)
end
label!(ax, x, y, s; align = (:left, :center), kw...) =
    text!(ax, x, y; text = s, align, fontsize = 14, color = C.ink, kw...)

# vertical guides
for (x, s) in ((-15, L"r_*/2"), (-14, L"r_*=2^{-14}"), (-12, L"a/2"), (-11, L"a=2^{-11}"))
    vlines!(axw, [x]; color = (C.ink2, 0.5), linestyle = :dot)
    vlines!(axr, [x]; color = (C.ink2, 0.5), linestyle = :dot)
end

y = 5.2   # Fock ball
bar!(axr, xr[1], -14, y, C.blue)
bar!(axr, -14, -13, y, C.blue; alpha = 0.3)
label!(
    axr,
    -12.8,
    y,
    L"Fock series $F_N$, all angles: tail $\le 2^{160}(N+2)^8 2^{-N}$ on $\rho\le r_*$ ($N=244$ for $2^{-20}$); pale: sharpened $r_*=2^{-13}$",
)
y = 4.2   # blend
bar!(axr, -15, -14, y, C.yellow)
label!(axr, -13.8, y, L"blend zone $[r_*/2,\,r_*]$: $\chi$ falls from 1 to 0 ($C^\infty$)")
y = 3.2   # extraction shell
bar!(axr, -12, -11, y, C.aqua)
label!(
    axr,
    -10.8,
    y,
    L"extraction shell $[a/2, a]$: trace inversions averaged over $r$ give the free data $h_n$",
)
y = 2.2   # global
bar!(axr, xr[1], -15, y, C.gray; alpha = 0.25)
bar!(axr, -15, xr[2], y, C.gray)
label!(
    axr,
    -9.5,
    y,
    L"global heat approximant $W$ (sup and $L^2$ certifiable everywhere; weight $1-\chi$)";
    color = :white,
)
y = 1.2   # clipping cores
for (k, M) in zip(epsk, Ms)
    scatter!(axr, [-l2(M)], [y]; color = C.violet, markersize = 12)
    text!(
        axr,
        -l2(M),
        y;
        text = L"\varepsilon=2^{-%$k}",
        align = (:center, :bottom),
        offset = (0, 8),
        fontsize = 12,
        color = C.ink2,
    )
end
label!(
    axr,
    -9.6,
    y,
    L"Coulomb clipping radius $1/M$ ($M = 36\|v\|/\varepsilon$; heat budget rows)",
)
y = 0.2   # physical scales
for (x, s) in (
    (l2(sqrt(T_HEAT)), L"\sqrt{t_*}=2^{-5}"),
    (l2(0.1), "blend radius of\nFig. patched_blend"),
    (l2(1.18), "empirical Fock\nscale ≈ 1.2 bohr"),
)
    scatter!(axr, [x], [y]; color = C.orange, marker = :diamond, markersize = 13)
    text!(
        axr,
        x,
        y;
        text = s,
        align = (:center, :top),
        offset = (0, -10),
        fontsize = 12,
        color = C.ink2,
    )
end
label!(axr, -21.5, y, "physical / heat scales →"; color = C.ink2)

# ------------------------------------------------------------------ pipeline
axp = Axis(
    fig[3, 1];
    height = 430,
    title = "(b) pipeline: green = executed and certified, orange dashed = specified with explicit error budget, never executed",
    titlesize = 17,
    titlefont = :regular,
)
hidedecorations!(axp)
hidespines!(axp)
limits!(axp, 0, 100, 0, 40)
function box!(x, y, w, h, s; done = false)
    col = done ? C.aqua : C.orange
    poly!(
        axp,
        Rect(x, y, w, h);
        color = (col, 0.12),
        strokecolor = col,
        strokewidth = 2,
        linestyle = done ? :solid : :dash,
    )
    text!(
        axp,
        x + w / 2,
        y + h / 2;
        text = s,
        align = (:center, :center),
        fontsize = 13,
        color = C.ink,
        justification = :center,
    )
end
arrow!(x0, y0, x1, y1) = arrows2d!(
    axp,
    [x0],
    [y0],
    [x1 - x0],
    [y1 - y0];
    color = C.ink2,
    shaftwidth = 1.5,
    tipwidth = 10,
    tiplength = 10,
)

res = parse(Float64, jsonstring(seed, "full_residual_norm_upper"))
lo = @sprintf(
    "%.10f",
    floor(parse(BigFloat, jsonstring(seed, "energy_lower")) * big(10)^10) / big(10)^10
)
hi = @sprintf(
    "%.10f",
    ceil(parse(BigFloat, jsonstring(seed, "energy_upper")) * big(10)^10) / big(10)^10
)
l2e = parse(Float64, jsonstring(seed, "oriented_normalized_L2_error_upper"))
box!(
    1,
    15,
    17,
    11,
    @sprintf("seed trial v\n151-term log-Hylleraas\nfull PDE residual\n≤ %.3g", res);
    done = true,
)
box!(
    22,
    15,
    18,
    11,
    @sprintf("Temple (β = −5/2)\nE₀ ∈ [%s,\n%s]\n‖v − ψ‖₂ ≤ %.2g", lo, hi, l2e);
    done = true,
)
arrow!(18, 20.5, 22, 20.5)
# global branch (top)
box!(45, 28, 16, 11, "heat smoothing\nexp(−t*H), t* = 1/1024\n‖·‖₂→∞ < 4096")
box!(65, 28, 17, 11, "clip 1/r at M, Dyson\norder J, Gaussian\ncubature (dim 7J+6)")
box!(86, 28, 13, 11, rich("global W", subscript("k"), "\nsup & L² error ≤ ε"))
arrow!(40, 23, 45, 31)
arrow!(61, 33.5, 65, 33.5)
arrow!(82, 33.5, 86, 33.5)
# local branch (bottom)
box!(
    45,
    2,
    16,
    12,
    rich(
        "averaged trace\ninversion, r ∈ [a/2, a]\nfree data hₙ: needs\nδ ≤ ε·2",
        superscript("−(56+12n)"),
    ),
)
box!(
    65,
    2,
    17,
    12,
    rich(
        "recurrence + C⁰ heat-\nintegral inverses → pₙ\n‖pₙ‖ ≤ 2",
        superscript("155"),
        " (n+1)",
        superscript("8"),
        " 8192",
        superscript("n"),
    ),
)
box!(
    86,
    2,
    13,
    12,
    rich("Fock F", subscript("N,k"), "\nN = 244 for\n2", superscript("−20"), " on ρ ≤ r*"),
)
arrow!(40, 18, 45, 10)
arrow!(61, 8, 65, 8)
arrow!(82, 8, 86, 8)
# blend
box!(
    84,
    16.5,
    15.5,
    9,
    rich(
        "Z",
        subscript("k"),
        " = χF",
        subscript("N,k"),
        " + (1−χ)W",
        subscript("k"),
        "\nglobal sup, L² error → 0",
    ),
)
arrow!(92.5, 28, 92.5, 25.5)
arrow!(92.5, 14, 92.5, 16.5)
text!(
    axp,
    1,
    4;
    text = "Budget arithmetic (EFFECTIVE_*_BUDGET.py) was executed and is re-derived\nin patched_budget.jl; the cubature and the coefficient evaluator were not.",
    align = (:left, :center),
    fontsize = 13,
    color = C.ink2,
)
rowgap!(fig.layout, 1, 4)
savefig("patched_schematic", fig)

ncgroup(
    "schematic";
    attrib = ["comment" => "Radii of the EFFECTIVE_* construction (bohr)"],
) do g
    putvar!(
        g,
        "radius_name_index",
        collect(1:6),
        ("radius",);
        attrib = ["long_name" => "1 r*/2, 2 r*, 3 r* sharpened, 4 a/2, 5 a, 6 sqrt(t*)"],
    )
    putvar!(
        g,
        "radius",
        Float64.([R_STAR / 2, R_STAR, R_STAR_SHARP, A_SHELL / 2, A_SHELL, sqrt(T_HEAT)]),
        ("radius",);
        attrib = ["units" => "bohr", "long_name" => "radius of the construction, named by index"],
    )
    putvar!(
        g,
        "clipping_M",
        Float64.(Ms),
        ("heat_row",);
        attrib = ["long_name" => "clipping level M"],
    )
    putvar!(
        g,
        "log2_inv_eps",
        Float64.(epsk),
        ("heat_row",);
        attrib = ["long_name" => "1/ε", view_hint(value_transform = "exp2")...],
    )
end
