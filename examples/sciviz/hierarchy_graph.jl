# Schematic of the Fock recurrence as a dependency graph on the (k, p) lattice:
#   (Λ² − k(k+4))ψ_kp = 2(k+2)(p+1)ψ_k,p+1 + (p+1)(p+2)ψ_k,p+2 − 2Vψ_k−1,p + 2Eψ_k−2,p.
# Run: julia --project=examples/sciviz examples/sciviz/hierarchy_graph.jl
include(joinpath(@__DIR__, "common.jl"))
const latexstring = CairoMakie.Makie.LaTeXStrings.latexstring

const KMAX = 6
const COL =
    (V = "#2a78d6", E = "#eb6834", log1 = "#1baf7a", log2 = "#e87ba4", solv = "#7a7a7a")

pos(k, p) = Point2f(k, -1.4p)

function edge!(ax, a, b; color, bend = 0.0, linestyle = :solid, linewidth = 2)
    pa, pb = pos(a...), pos(b...)
    d = pb - pa
    n = Point2f(-d[2], d[1])
    mid = (pa + pb) / 2 + bend * n
    ts = range(0, 1; length = 40)
    curve = [(1 - t)^2 * pa + 2t * (1 - t) * mid + t^2 * pb for t in ts]
    # shorten so arrows stop at the node circle
    r = 0.2f0
    keep = [norm2(c - pa) > r && norm2(c - pb) > r for c in curve]
    pts = curve[keep]
    lines!(ax, pts; color, linestyle, linewidth)
    tip = pts[end]
    dir = tip - pts[end-1]
    scatter!(
        ax,
        [tip];
        marker = :rtriangle,
        markersize = 11,
        color,
        rotation = atan(dir[2], dir[1]),
    )
end
norm2(v) = sqrt(sum(abs2, v))

fig = Figure(size = (1100, 640), fontsize = 16)
ax = Axis(
    fig[1, 1];
    aspect = DataAspect(),
    xticks = 0:KMAX,
    xlabel = L"k\ \ (\mathrm{power\ of\ } R)",
    ylabel = L"p\ \ (\mathrm{power\ of\ } \log R)",
    yticks = (-1.4 .* (0:(KMAX÷2)), [latexstring("$p") for p = 0:(KMAX÷2)]),
    xgridvisible = false,
    ygridvisible = false,
    rightspinevisible = false,
    topspinevisible = false,
)
nodes = [(k, p) for k = 0:KMAX for p = 0:(k÷2)]
# even k is resonant: Λ² − k(k+4) annihilates the harmonics Y_kl
for k = 2:2:KMAX
    poly!(
        ax,
        Rect2f(k - 0.42, -1.4 * (KMAX ÷ 2) - 0.75, 0.84, 1.4 * (KMAX ÷ 2) + 1.5);
        color = (:gray, 0.12),
        strokewidth = 0,
    )
    text!(
        ax,
        Point2f(k, -1.4 * (KMAX ÷ 2) - 0.7);
        text = "resonant",
        align = (:center, :bottom),
        fontsize = 12,
        color = :gray35,
    )
end
for (k, p) in nodes
    p <= (k - 1) ÷ 2 && k >= 1 && edge!(ax, (k - 1, p), (k, p); color = COL.V)
    p <= (k - 2) ÷ 2 && k >= 2 && edge!(ax, (k - 2, p), (k, p); color = COL.E, bend = 0.25)
    p + 1 <= k ÷ 2 && edge!(ax, (k, p + 1), (k, p); color = COL.log1, bend = -0.18)
    p + 2 <= k ÷ 2 && edge!(ax, (k, p + 2), (k, p); color = COL.log2, bend = -0.32)
    # even k: solvability of the p equation fixes the Y_kl part of ψ_k,p+1
    iseven(k) &&
        k >= 2 &&
        p + 1 <= k ÷ 2 &&
        edge!(
            ax,
            (k, p),
            (k, p + 1);
            color = COL.solv,
            bend = -0.18,
            linestyle = :dash,
            linewidth = 1.5,
        )
end
for (k, p) in nodes
    free = iseven(k) && p == 0 && k >= 2
    scatter!(
        ax,
        [pos(k, p)];
        markersize = 46,
        color = free ? "#fff4d6" : :white,
        strokecolor = :black,
        strokewidth = 1.5,
    )
    text!(
        ax,
        pos(k, p);
        text = latexstring("\\psi_{$k$p}"),
        align = (:center, :center),
        fontsize = 17,
    )
end
freelabels = Dict(2 => L"a_{21}", 4 => L"a_{40},a_{42}", 6 => L"a_{6l}")
for (k, s) in freelabels
    text!(
        ax,
        pos(k, 0) + Point2f(0, 0.42);
        text = s,
        align = (:center, :bottom),
        fontsize = 14,
        color = "#8a5a00",
    )
end
xlims!(ax, -0.6, KMAX + 0.6)
ylims!(ax, -1.4 * (KMAX ÷ 2) - 0.9, 1.0)
elems = [
    LineElement(color = COL.V, linewidth = 2),
    LineElement(color = COL.E, linewidth = 2),
    LineElement(color = COL.log1, linewidth = 2),
    LineElement(color = COL.log2, linewidth = 2),
    LineElement(color = COL.solv, linewidth = 1.5, linestyle = :dash),
    MarkerElement(
        marker = :circle,
        color = "#fff4d6",
        strokecolor = :black,
        strokewidth = 1.5,
        markersize = 18,
    ),
    PolyElement(color = (:gray, 0.12)),
]
labels = [
    L"-2V\psi_{k-1,p}\ \ \mathrm{(Coulomb)}",
    L"+2E\psi_{k-2,p}\ \ \mathrm{(energy)}",
    L"2(k+2)(p+1)\psi_{k,p+1}",
    L"(p+1)(p+2)\psi_{k,p+2}",
    L"\mathrm{solvability\ fixes\ } Y_{kl}\mathrm{\ part\ (even\ } k)",
    L"\mathrm{free\ harmonic\ constants\ in\ }\psi_{k,0}",
    L"\mathrm{resonant\ } k: \ (\Lambda^2-k(k+4))Y_{kl}=0",
]
Legend(
    fig[1, 2],
    elems,
    labels,
    L"\mathrm{source\ of\ }(\Lambda^2-k(k+4))\psi_{k,p}";
    framevisible = false,
    titlesize = 15,
    labelsize = 14,
)
savefig("hierarchy_graph", fig)
