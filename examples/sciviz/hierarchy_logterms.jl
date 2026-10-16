# Logarithmic coefficients of the Fock expansion.
#  hierarchy_logterms:  ψ21, ψ31, ψ41, ψ42 from FockExpansion.jl; ψ41 split by powers of Z
#    (Zψ⁽¹⁾ and Z³ψ⁽³⁾ closed forms, Z²ψ⁽²⁾ by quadrature, and the harmonic part c₀Y₄₀ + c₂Y₄₂
#    fixed by solvability of the ψ40 equation); and c₀, c₂ as functions of a₂₁ and E.
#  hierarchy_leadinglog: the leading-log tower ψ_{2m,m} from the spectral chain. Each is a
#    single degree-2m harmonic and does not depend on E or on the free constants.
# Needs output/hierarchy_package.nc from hierarchy_compute.jl.
# Run: julia -t 4 --project=examples/sciviz examples/sciviz/hierarchy_logterms.jl
include(joinpath(@__DIR__, "hierarchy_common.jl"))

const FE = FockExpansion
pk = NCDataset(joinpath(OUT, "hierarchy_package.nc")) do ds
    (
        α = Array(ds["alpha"]),
        θ = Array(ds["theta"]),
        f = Dict(name => Array(ds[name]) for name in ("psi21", "psi31", "psi41", "psi42")),
    )
end
αs, θs = pk.α, pk.θ
Z, E, a21 = HE.Z, HE.E, HE.a21

# ψ41 by powers of Z (src/fourth_order.jl): the Z² part is ψ41 minus the closed-form parts.
ξ(a, t) = sqrt(1 - sin(a) * cos(t))
z1(a, t) = Z * (π - 2) / (2880π) * (3 * (32E - 15) - 8 * (12E - 5) * ξ(a, t)^2)
z3(a, t) = -Z^3 * (π - 2) / (120π) * (4 + 5sin(a)) * sin(a) * cos(t)
harm(a, t) = FE.psi41_harmonic(a, t; Z, E, a21)
P1 = grid(z1, αs, θs)
P3 = grid(z3, αs, θs)
PH = grid(harm, αs, θs)
P2 = pk.f["psi41"] .- P1 .- P3 .- PH

ticks = ([0, π / 2, π], [L"0", L"\pi/2", L"\pi"])
function panel!(pos, A, title; first = false, xlab = true)
    m = maximum(abs, A)
    ax = Axis(
        pos;
        aspect = DataAspect(),
        xticks = ticks,
        yticks = ticks,
        title = latexstring(title, "\\quad", @sprintf("%.2g", m)),
        yticklabelsvisible = first,
        xticklabelsvisible = xlab,
        xlabel = xlab ? L"\alpha" : "",
        ylabel = first ? L"\theta" : "",
        xticklabelsize = 12,
        yticklabelsize = 12,
        titlesize = 15,
    )
    heatmap!(ax, αs, θs, A ./ m; colormap = DIVERGING, colorrange = (-1, 1))
    minimum(A) < 0 < maximum(A) &&
        contour!(ax, αs, θs, A; levels = [0.0], color = :black, linewidth = 0.8)
    ax
end

fig = Figure(size = (1400, 680), fontsize = 15)
panel!(fig[1, 1], pk.f["psi21"], "\\psi_{2,1}"; first = true, xlab = false)
panel!(fig[1, 2], pk.f["psi31"], "\\psi_{3,1}"; xlab = false)
panel!(fig[1, 3], pk.f["psi41"], "\\psi_{4,1}"; xlab = false)
panel!(fig[1, 4], pk.f["psi42"], "\\psi_{4,2}"; xlab = false)
panel!(fig[2, 1], P1, "Z\\,\\psi^{(1)}_{4,1}"; first = true)
panel!(fig[2, 2], P2, "Z^2\\psi^{(2)}_{4,1}")
panel!(fig[2, 3], P3, "Z^3\\psi^{(3)}_{4,1}")
panel!(fig[2, 4], PH, "c_0Y_{40}+c_2Y_{42}")
Label(
    fig[1, 0],
    L"\mathrm{log\ terms}";
    rotation = π / 2,
    tellheight = false,
    fontsize = 16,
)
Label(
    fig[2, 0],
    L"\psi_{4,1}\ \mathrm{by\ powers\ of\ }Z";
    rotation = π / 2,
    tellheight = false,
    fontsize = 16,
)
Colorbar(
    fig[1:2, 5];
    colormap = DIVERGING,
    colorrange = (-1, 1),
    label = L"\psi/\max|\psi|\ \ (\max\ \mathrm{in\ title})",
    height = Relative(0.7),
)

# Harmonic part versus a21 and E: c = Z[(3π−8)E/(540π) − (5π−14)a21/(90π)] adds c(Y40 + 4Y42).
c0(Zv, Ev, av) =
    evalpoly(Zv, (0, FE.PSI41_HARMONIC.l0...)) +
    Zv * ((3π - 8) * Ev / (540π) - (5π - 14) * av / (90π))
c2(Zv, Ev, av) =
    evalpoly(Zv, (0, FE.PSI41_HARMONIC.l2...)) +
    4Zv * ((3π - 8) * Ev / (540π) - (5π - 14) * av / (90π))
ax = Axis(
    fig[1:2, 6];
    xlabel = L"a_{21}",
    ylabel = L"\mathrm{coefficient\ in\ }\psi_{4,1}",
    title = L"\mathrm{harmonic\ part\ of\ }\psi_{4,1}\ (Z=2)",
)
as = range(0, 1; length = 50)
for (i, Ev) in enumerate((-2.0, E, -4.0))
    lab = i == 2 ? L"E_{\mathrm{He}}" : latexstring(@sprintf("E=%.0f", Ev))
    lines!(
        ax,
        as,
        c0.(Z, Ev, as);
        color = SERIES[i],
        linewidth = 2,
        label = latexstring("c_0,\\ ", lab),
    )
    lines!(
        ax,
        as,
        c2.(Z, Ev, as);
        color = SERIES[i],
        linewidth = 2,
        linestyle = :dash,
        label = latexstring("c_2,\\ ", lab),
    )
end
scatter!(ax, [a21, a21], [c0(Z, E, a21), c2(Z, E, a21)]; color = :black, markersize = 10)
text!(
    ax,
    a21,
    c0(Z, E, a21);
    text = L"\ \ \mathrm{He}",
    align = (:left, :top),
    fontsize = 13,
)
text!(
    ax,
    a21,
    c2(Z, E, a21);
    text = L"\ \ \mathrm{He}",
    align = (:left, :bottom),
    fontsize = 13,
)
axislegend(ax; position = :lb, framevisible = false, labelsize = 12, nbanks = 2)
colsize!(fig.layout, 6, Relative(0.26))
savefig("hierarchy_logterms", fig)
@printf(
    "He: c0 = %.6f, c2 = %.6f; max|Z²ψ⁽²⁾| = %.3e\n",
    c0(Z, E, a21),
    c2(Z, E, a21),
    maximum(abs, P2)
)

# --- leading-log tower ψ_{2m,m} ---
const MMAX = 6
g, ψ = he_chain(; n = 64, kmax = 2MMAX)
gE, ψE = let gg = SGrid(Float64, 64)
    gg,
    sfock(
        gg;
        Z,
        E = -1.0,
        kmax = 2MMAX,
        free = Dict((2, 1) => 0.1, (4, 0) => 0.3, (4, 2) => -0.2, (6, 1) => 0.05),
    )
end
tower = [(2m, m) for m = 1:MMAX]
lα = interior(121)
vals = interp_grid(g, ψ, tower, lα, lα)
fig = Figure(size = (1400, 600), fontsize = 15)
top = fig[1, 1] = GridLayout()
bot = fig[2, 1] = GridLayout()
for (i, key) in enumerate(tower)
    A = vals[key]
    m = maximum(abs, A)
    local ax = Axis(
        top[1, i];
        aspect = DataAspect(),
        xticks = ticks,
        yticks = ticks,
        yticklabelsvisible = i == 1,
        title = latexstring("\\psi_{$(key[1]),$(key[2])}"),
        xlabel = L"\alpha",
        ylabel = i == 1 ? L"\theta" : "",
    )
    heatmap!(ax, lα, lα, A ./ m; colormap = DIVERGING, colorrange = (-1, 1))
    contour!(ax, lα, lα, A; levels = [0.0], color = :black, linewidth = 0.8)
end
Colorbar(
    top[1, MMAX+1];
    colormap = DIVERGING,
    colorrange = (-1, 1),
    label = L"\psi/\max|\psi|",
    height = Relative(0.6),
)
# size, harmonic purity and independence of E and the free constants
ms = 1:MMAX
size_ = [rms(g, ψ[(2m, m)]) for m in ms]
purity = Float64[]
for m in ms
    f = ψ[(2m, m)]
    Hs = [Y ./ rms(g, Y) for (l, Y) in singlet_harmonics(g, 2m)]
    W = g.W ./ sum(g.W)
    proj = sum(sum(W .* f .* Y) * Y for Y in Hs)
    push!(purity, rms(g, f - proj) / rms(g, f))
end
indep = [rms(g, ψ[(2m, m)] - ψE[(2m, m)]) / rms(g, ψ[(2m, m)]) for m in ms]
for m in ms
    @printf(
        "ψ%d%d: rms %.3e, non-harmonic fraction %.1e, change with E and free data %.1e\n",
        2m,
        m,
        size_[m],
        purity[m],
        indep[m]
    )
end
ax = Axis(
    bot[1, 1];
    yscale = log10,
    xlabel = L"m",
    ylabel = L"\mathrm{RMS\ of\ }\psi_{2m,m}",
    xticks = ms,
    yticks = LogTicks(-14:3:0),
    title = L"\mathrm{(b)\ size\ of\ the\ leading\ log}",
)
scatterlines!(ax, ms, size_; color = SERIES[1], linewidth = 2, markersize = 9)
ax2 = Axis(
    bot[1, 2];
    yscale = log10,
    xlabel = L"m",
    xticks = ms,
    ylabel = L"\mathrm{relative\ size}",
    yticks = LogTicks(-17:2:-8),
    title = L"\mathrm{(c)\ a\ single\ harmonic,\ independent\ of\ }E\mathrm{\ and\ free\ data}",
)
scatterlines!(
    ax2,
    ms,
    max.(purity, 1e-17);
    color = SERIES[2],
    linewidth = 2,
    markersize = 9,
    label = L"\mathrm{part\ outside\ degree\ }2m\mathrm{\ harmonics}",
)
scatterlines!(
    ax2,
    ms,
    max.(indep, 1e-17);
    color = SERIES[3],
    linewidth = 2,
    markersize = 9,
    label = L"\mathrm{change\ for\ }E=-1,\ a_{21}=0.1,\ a_{4l},a_{6l}\neq0",
)
ylims!(ax2, 1e-17, 1e-8)
axislegend(ax2; position = :lt, framevisible = false, labelsize = 13)
rowsize!(fig.layout, 2, Relative(0.42))
Label(
    top[0, 1:MMAX],
    L"\mathrm{(a)\ leading\ log\ tower\ }\psi_{2m,m}\mathrm{\ (coefficient\ of\ }R^{2m}\log^m R)";
    fontsize = 15,
)
savefig("hierarchy_leadinglog", fig)
