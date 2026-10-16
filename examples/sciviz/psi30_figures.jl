# Paper figures for ψ30, the third-order nonlogarithmic Fock coefficient.
#
#     julia -t 4 --project=examples/sciviz examples/sciviz/psi30_compute.jl   # once, caches grids
#     julia -t 4 --project=examples/sciviz examples/sciviz/psi30_figures.jl
#
# Reads output/psi30.nc and output/psi30_agreement.nc; one-dimensional data and the
# zoomed coalescence grid are cheap and computed here.
include("common.jl")
using ForwardDiff, QuadGK, Printf

# ---------------------------------------------------------------------------- loading

function load_angular(path)
    NCDataset(path) do ds
        αs = Float64.(ds["alpha"][:])
        θs = Float64.(ds["theta"][:])
        names = [k for k in keys(ds) if !(k in ("alpha", "theta", "x", "y", "z"))]
        fields = Dict(k => Float64.(ds[k][:, :]) for k in names)
        (; αs, θs, fields)
    end
end

const MAIN = load_angular(joinpath(OUT, "psi30.nc"))
const AGREE = load_angular(joinpath(OUT, "psi30_agreement.nc"))
const αs, θs = MAIN.αs, MAIN.θs
F(name) = MAIN.fields[name]

# ---------------------------------------------------------------------------- helpers

const PITICKS = (
    [0, π / 4, π / 2, 3π / 4, π],
    [L"0", L"\pi/4", L"\pi/2", L"3\pi/4", L"\pi"],
)
const DIVERGING = Reverse(:RdBu)

function angle_axis(pos; kw...)
    Axis(
        pos;
        xlabel = L"\alpha",
        ylabel = L"\theta",
        xticks = PITICKS,
        yticks = PITICKS,
        aspect = DataAspect(),
        limits = (0, π, 0, π),
        kw...,
    )
end

"""Diverging heatmap centred on zero, with optional level lines and nodal line."""
function signed_map!(ax, xs, ys, A; levels = 12, nodal = true, colorscale = identity)
    lo, hi = extrema(A)
    signed = lo < 0 < hi
    m = maximum(abs, A)
    # Diverging colours only for fields that change sign; one-signed fields get viridis.
    hm = heatmap!(
        ax,
        xs,
        ys,
        A;
        colormap = signed ? DIVERGING : :viridis,
        colorrange = signed ? (-m, m) : (lo, hi),
        colorscale,
        rasterize = 2,
    )
    if levels !== nothing
        contour!(ax, xs, ys, A; levels, color = (:black, 0.35), linewidth = 0.6)
    end
    if nodal && minimum(A) < 0 < maximum(A)
        contour!(ax, xs, ys, A; levels = [0.0], color = :black, linewidth = 1.6)
    end
    hm
end

"""Sequential log-colour heatmap of a non-negative field."""
function log_map!(ax, xs, ys, A; floor = 1e-17, colormap = :magma, colorrange = nothing)
    B = max.(A, floor)
    cr = colorrange === nothing ? (max(minimum(B), floor), maximum(B)) : colorrange
    heatmap!(ax, xs, ys, B; colormap, colorscale = log10, colorrange = cr, highclip = :white, rasterize = 2)
end

"""Mark the two-particle coalescences: r₂ = 0 on α = 0, r₁ = 0 on α = π
(r₁/R = cos(α/2), r₂/R = sin(α/2)), and r₁₂ = 0 at (α, θ) = (π/2, 0)."""
function mark_coalescences!(ax; labels = true)
    vlines!(ax, [0.0, π]; color = :seagreen, linewidth = 5)
    scatter!(
        ax,
        [π / 2],
        [0.0];
        marker = :star5,
        markersize = 22,
        color = :gold,
        strokecolor = :black,
        strokewidth = 1.2,
    )
    if labels
        text!(ax, π / 2, 0.12; text = L"r_{12}=0", align = (:center, :bottom), fontsize = 15)
        text!(ax, 0.1, π / 2; text = L"r_2=0", rotation = π / 2, align = (:center, :top), fontsize = 15)
        text!(ax, π - 0.1, π / 2; text = L"r_1=0", rotation = π / 2, align = (:center, :bottom), fontsize = 15)
    end
end

"""Hyperspherical volume weight sin²α sinθ dα dθ on the cached grid."""
const WEIGHT = [sin(a)^2 * sin(t) for a in αs, t in θs] .* (step(range(αs[1], αs[end]; length = length(αs)))^2)
wnorm(A) = sqrt(sum(WEIGHT .* A .^ 2) / sum(WEIGHT))

"""Shrink the figure to its layout (removes the slack left by DataAspect axes), then save."""
function save_tight(name, fig)
    resize_to_layout!(fig)
    savefig(name, fig)
end

# ---------------------------------------------------------------------------- figures

function fig_map()
    ψ = F("psi30")
    fig = Figure(size = (720, 620))
    ax = angle_axis(fig[1, 1]; title = L"\psi_{3,0}(\alpha,\theta)\ \ \mathrm{(helium,\ ground\ state)}")
    hm = signed_map!(ax, αs, θs, ψ; levels = 16)
    mark_coalescences!(ax)
    Colorbar(fig[1, 2], hm; label = L"\psi_{3,0}")
    save_tight("psi30_map", fig)
end

function fig_hypersphere()
    ψ = F("psi30")
    X = [hypersphere_xyz(a, t)[1] for a in αs, t in θs]
    Y = [hypersphere_xyz(a, t)[2] for a in αs, t in θs]
    Zc = [hypersphere_xyz(a, t)[3] for a in αs, t in θs]
    cr = extrema(ψ)
    fig = Figure(size = (1200, 560))
    ax1 = Axis3(
        fig[1, 1];
        aspect = :data,
        azimuth = -0.2π,
        elevation = 0.22π,
        xlabel = L"r_1/R",
        ylabel = L"(r_2/R)\cos\theta",
        zlabel = L"(r_2/R)\sin\theta",
        title = L"\psi_{3,0}\ \mathrm{on\ the\ shape\ hemisphere}",
    )
    sf = surface!(ax1, X, Y, Zc; color = ψ, colormap = :viridis, colorrange = cr, shading = NoShading, rasterize = 2)
    # Coalescence loci: r₂ = 0 is the pole (1,0,0); r₁ = 0 is the arc x = 0; r₁₂ = 0 is a point.
    scatter!(ax1, [1.0], [0.0], [0.0]; color = :seagreen, markersize = 18)
    ts = range(0, π; length = 100)
    lines!(ax1, zeros(100), cos.(ts), sin.(ts); color = :seagreen, linewidth = 4)
    scatter!(ax1, [cos(π / 4)], [sin(π / 4)], [0.0]; marker = :star5, markersize = 24,
        color = :gold, strokecolor = :black, strokewidth = 1.2)
    text!(ax1, Point3f(1.0, 0.0, 0.08); text = L"r_2=0", fontsize = 15)
    text!(ax1, Point3f(0.0, -0.1, 1.05); text = L"r_1=0", fontsize = 15)
    text!(ax1, Point3f(cos(π / 4) + 0.05, sin(π / 4) - 0.05, 0.12); text = L"r_{12}=0", fontsize = 15)

    ax2 = Axis3(
        fig[1, 2];
        azimuth = 1.3π,
        elevation = 0.2π,
        xlabel = L"\alpha",
        ylabel = L"\theta",
        zlabel = L"\psi_{3,0}",
        xticks = PITICKS,
        yticks = PITICKS,
        title = L"\mathrm{height\ field\ over\ }(\alpha,\theta)",
    )
    surface!(ax2, αs, θs, ψ; color = ψ, colormap = :viridis, colorrange = cr, rasterize = 2)
    Colorbar(fig[1, 3], sf; label = L"\psi_{3,0}", height = Relative(0.7))
    colgap!(fig.layout, 1, 90)
    savefig("psi30_hypersphere", fig)
end

function fig_parts()
    fig = Figure(size = (1100, 900))
    names = [
        ("psi30_classical", L"\mathrm{(a)\ dilogarithmic\ part\ }\mathcal{D}"),
        ("psi30_state", L"\mathrm{(b)\ }E,\,a_{21}\mathrm{\ part}"),
        ("psi30_elliptic", L"\mathrm{(c)\ elliptic\ part\ }\sum F"),
    ]
    for (k, (name, title)) in enumerate(names)
        r, c = fldmod1(k, 2)
        ax = angle_axis(fig[r, 2c-1]; title)
        hm = signed_map!(ax, αs, θs, F(name); levels = 12)
        Colorbar(fig[r, 2c], hm)
    end
    ax = angle_axis(fig[2, 3]; title = L"\mathrm{(d)\ }\log_{10}\!\left(\sum|\mathrm{parts}|/|\psi_{3,0}|\right)")
    C = F("log10_cancellation")
    hm = heatmap!(ax, αs, θs, C; colormap = :inferno, colorrange = (0, maximum(C)), rasterize = 2)
    contour!(ax, αs, θs, C; levels = 0.05:0.05:0.4, color = (:white, 0.6), linewidth = 0.6)
    Colorbar(fig[2, 4], hm; label = L"\mathrm{digits\ lost}")
    save_tight("psi30_parts", fig)
    for name in ("psi30_classical", "psi30_state", "psi30_elliptic", "psi30")
        @printf("%-18s weighted rms = %.4e   max|.| = %.4e\n", name, wnorm(F(name)), maximum(abs, F(name)))
    end
end

function fig_sectors()
    fig = Figure(size = (1100, 900))
    for n = 0:3
        r, c = fldmod1(n + 1, 2)
        title = L"\mathrm{coefficient\ of\ } Z^{%$n}"
        ax = angle_axis(fig[r, 2c-1]; title)
        hm = signed_map!(ax, αs, θs, F("sector_Z$n"); levels = 12)
        Colorbar(fig[r, 2c], hm)
    end
    save_tight("psi30_Z_sectors", fig)

    # Share of each Z^n sector in the hyperspherical L² norm as a function of Z.
    norms = [wnorm(F("sector_Z$n")) for n = 0:3]
    Zs = range(0.05, 6; length = 300)
    shares = [[z^n * norms[n+1] for z in Zs] for n = 0:3]
    total = sum(shares)
    fig = Figure(size = (720, 440))
    ax = Axis(fig[1, 1]; xlabel = L"Z", ylabel = L"\mathrm{share\ of\ the\ norm}",
        title = L"\mathrm{weight\ of\ the\ }Z^n\mathrm{\ sectors\ of\ }\psi_{3,0}\ (E=a_{21}=0)")
    colors = Makie.wong_colors()
    for n = 0:3
        lines!(ax, Zs, shares[n+1] ./ total; color = colors[n+1], linewidth = 2.5, label = L"n=%$n")
    end
    for (z, lab) in ((1, L"\mathrm{H}^-"), (2, L"\mathrm{He}"), (3, L"\mathrm{Li}^+"))
        vlines!(ax, [z]; color = (:gray, 0.6), linestyle = :dash)
        text!(ax, z, 1.0; text = lab, align = (:center, :bottom), fontsize = 14)
    end
    ylims!(ax, 0, 1.1)
    axislegend(ax; position = :rc)
    savefig("psi30_Z_shares", fig)
    @printf("sector norms ‖A_n‖ = %s\n", string(norms))
end

function fig_sensitivity()
    fig = Figure(size = (1500, 470))
    names = [
        ("dpsi30_dE", L"\mathrm{(a)\ }\partial\psi_{3,0}/\partial E"),
        ("dpsi30_da21", L"\mathrm{(b)\ }\partial\psi_{3,0}/\partial a_{21}"),
        ("dpsi30_dZ", L"\mathrm{(c)\ }\partial\psi_{3,0}/\partial Z\ \ (E,a_{21}\mathrm{\ fixed})"),
    ]
    for (k, (name, title)) in enumerate(names)
        ax = angle_axis(fig[1, 2k-1]; title)
        hm = signed_map!(ax, αs, θs, F(name); levels = 12)
        Colorbar(fig[1, 2k], hm)
    end
    save_tight("psi30_sensitivity", fig)
end

function fig_slices()
    Z, E, a21 = HE
    ε = 1e-4
    as = range(ε, π - ε; length = 1201)
    ts = range(ε, π - ε; length = 1201)
    rows = [
        (L"\mathrm{(a)\ collinear,\ }\theta\to0", L"\alpha", [(a, ε) for a in as], as),
        (L"\mathrm{(b)\ collinear,\ }\theta\to\pi", L"\alpha", [(a, π - ε) for a in as], as),
        (L"\mathrm{(c)\ }r_1=r_2\ (\alpha=\pi/2)", L"\theta", [(π / 2, t) for t in ts], ts),
    ]
    fig = Figure(size = (1500, 470))
    colors = Makie.wong_colors()
    for (k, (title, xl, pts, xs)) in enumerate(rows)
        P = [psi30_parts(p...; Z, E, a21) for p in pts]
        ax = Axis(fig[1, k]; title, xlabel = xl, xticks = PITICKS)
        hlines!(ax, [0.0]; color = (:gray, 0.5))
        lines!(ax, xs, getfield.(P, :classical); color = colors[1], linestyle = :dash, label = L"\mathcal{D}")
        lines!(ax, xs, getfield.(P, :state); color = colors[2], linestyle = :dash, label = L"E,a_{21}")
        lines!(ax, xs, getfield.(P, :elliptic); color = colors[3], linestyle = :dash, label = L"\sum F")
        lines!(ax, xs, getfield.(P, :value); color = :black, linewidth = 2.5, label = L"\psi_{3,0}")
    end
    Legend(fig[2, 1:3], fig.content[1]; orientation = :horizontal, framevisible = false)
    savefig("psi30_slices", fig)
end

"""Kato cusp checks. e–e: along α = π/2, dψ30/dξ → ψ20(coalescence)/2 with ξ = r12/R.
e–nucleus: the θ-average of ψ30 satisfies d⟨ψ30⟩/ds → −Z ψ20(r2 = 0) with s = r2/R."""
function fig_cusps()
    Z, E, a21 = HE
    ψ20c = psi20(π / 2, 1e-10; HE...)
    ψ20n = psi20(1e-10, 1.0; HE...)
    katoee = ψ20c / 2
    katoen = -Z * ψ20n

    ξs = 10 .^ range(-6, log10(0.9); length = 120)
    θof(ξ) = acos(1 - ξ^2)                       # α = π/2 ⇒ ξ² = 1 − cos θ
    dee = map(ξs) do ξ
        θ = 2asin(ξ / sqrt(2))
        dψdθ = ForwardDiff.derivative(t -> psi30(π / 2, t; HE...), θ)
        dψdθ * 2ξ / sin(θ)
    end

    nodes, weights = gauss(64, 0, π)
    avg(s) = sum(w * psi30(2asin(s), t; HE...) * sin(t) / 2 for (t, w) in zip(nodes, weights))
    ss = 10 .^ range(-4, log10(0.7); length = 80)
    den = [ForwardDiff.derivative(avg, s) for s in ss]

    fig = Figure(size = (1100, 820))
    ax1 = Axis(fig[1, 1]; xscale = log10, xlabel = L"\xi=r_{12}/R\ \ (\alpha=\pi/2)",
        ylabel = L"\partial\psi_{3,0}/\partial\xi", title = L"\mathrm{(a)\ electron\!-\!electron\ cusp}")
    lines!(ax1, ξs, dee; color = :black, linewidth = 2)
    hlines!(ax1, [katoee]; color = :firebrick, linestyle = :dash, label = L"\psi_{2,0}(\mathrm{coal.})/2")
    axislegend(ax1; position = :lb)
    ax2 = Axis(fig[2, 1]; xscale = log10, yscale = log10, xlabel = L"\xi",
        ylabel = L"|\partial_\xi\psi_{3,0}-\psi_{2,0}/2|")
    lines!(ax2, ξs, abs.(dee .- katoee); color = :black, linewidth = 2)
    ref = ξs .* abs(dee[60] - katoee) / ξs[60]
    lines!(ax2, ξs, ref; color = :gray, linestyle = :dot, label = L"\propto\xi")
    axislegend(ax2; position = :rb)

    ax3 = Axis(fig[1, 2]; xscale = log10, xlabel = L"s=r_2/R",
        ylabel = L"\partial\langle\psi_{3,0}\rangle_\theta/\partial s",
        title = L"\mathrm{(b)\ electron\!-\!nucleus\ cusp}")
    lines!(ax3, ss, den; color = :black, linewidth = 2)
    hlines!(ax3, [katoen]; color = :firebrick, linestyle = :dash, label = L"-Z\,\psi_{2,0}(r_2=0)")
    axislegend(ax3; position = :lb)
    ax4 = Axis(fig[2, 2]; xscale = log10, yscale = log10, xlabel = L"s",
        ylabel = L"|\partial_s\langle\psi_{3,0}\rangle+Z\psi_{2,0}|")
    lines!(ax4, ss, abs.(den .- katoen); color = :black, linewidth = 2)
    ref = ss .* abs(den[40] - katoen) / ss[40]
    lines!(ax4, ss, ref; color = :gray, linestyle = :dot, label = L"\propto s")
    axislegend(ax4; position = :rb)
    savefig("psi30_cusps", fig)
    @printf("Kato e-e: target %.10f, smallest-ξ value %.10f\n", katoee, dee[1])
    @printf("Kato e-n: target %.10f, smallest-s value %.10f\n", katoen, den[1])
end

"""Zoom at r12 → 0: ψ30 has a cone (|ξ| cusp); subtracting c0 + ψ20(coal.)ξ/2 leaves a
much smoother remainder."""
function fig_zoom()
    h = 0.15
    as = range(π / 2 - h, π / 2 + h; length = 241)
    ts = range(1e-4, 2h; length = 241)
    Ψ = grid((a, t) -> psi30(a, t; HE...), as, ts)
    ξ = [sqrt(1 - sin(a) * cos(t)) for a in as, t in ts]
    ψ20c = psi20(π / 2, 1e-10; HE...)
    t0 = 1e-6
    c0 = psi30(π / 2, t0; HE...) - ψ20c / 2 * sqrt(1 - cos(t0))
    rem = Ψ .- c0 .- ψ20c / 2 .* ξ
    @printf("coalescence value ψ30(π/2, 0) ≈ %.12f\n", c0)

    fig = Figure(size = (1250, 470))
    zoomaxis(pos, title) = Axis(pos; title, xlabel = L"\alpha-\pi/2", ylabel = L"\theta", aspect = DataAspect())
    xs = as .- π / 2
    ax1 = zoomaxis(fig[1, 1], L"\mathrm{(a)\ }\psi_{3,0}-\psi_{3,0}(\mathrm{coal.})")
    hm1 = signed_map!(ax1, xs, ts, Ψ .- c0; levels = 14, nodal = false)
    scatter!(ax1, [0.0], [0.0]; marker = :star5, markersize = 18, color = :gold, strokecolor = :black, strokewidth = 1)
    Colorbar(fig[1, 2], hm1)
    ax2 = zoomaxis(fig[1, 3], L"\mathrm{(b)\ minus\ Kato\ cone\ }\frac{1}{2}\psi_{2,0}\,\xi")
    hm2 = signed_map!(ax2, xs, ts, rem; levels = 14, nodal = true)
    scatter!(ax2, [0.0], [0.0]; marker = :star5, markersize = 18, color = :gold, strokecolor = :black, strokewidth = 1)
    Colorbar(fig[1, 4], hm2)
    ax3 = zoomaxis(fig[1, 5], L"\mathrm{(c)\ }|\mathrm{remainder}|/\xi^2")
    hm3 = heatmap!(ax3, xs, ts, abs.(rem) ./ ξ .^ 2; colormap = :viridis, rasterize = 2)
    Colorbar(fig[1, 6], hm3)
    save_tight("psi30_coalescence_zoom", fig)
end

function fig_agreement()
    a, t = AGREE.αs, AGREE.θs
    dg = AGREE.fields["abs_diff_green"]
    dl = AGREE.fields["abs_diff_langner"]
    cr = (1e-16, max(maximum(dg), maximum(dl)))
    fig = Figure(size = (1150, 500))
    ax1 = angle_axis(fig[1, 1]; title = L"\mathrm{(a)\ }|\psi_{3,0}-\psi_{3,0}^{\mathrm{Green}}|")
    hm = log_map!(ax1, a, t, dg; colorrange = cr, floor = 1e-16)
    ax2 = angle_axis(fig[1, 2]; title = L"\mathrm{(b)\ }|\psi_{3,0}-\psi_{3,0}^{\mathrm{Langner}}|")
    log_map!(ax2, a, t, dl; colorrange = cr, floor = 1e-16)
    Colorbar(fig[1, 3], hm; label = L"\mathrm{absolute\ difference}")
    save_tight("psi30_agreement", fig)
    @printf("Green:   median %.2e  max %.2e\n", sort(vec(dg))[end÷2], maximum(dg))
    @printf("Langner: median %.2e  max %.2e\n", sort(vec(dl))[end÷2], maximum(dl))
end

function fig_derivatives()
    fig = Figure(size = (1150, 500))
    ax1 = angle_axis(fig[1, 1]; title = L"\mathrm{(a)\ }|\nabla\psi_{3,0}|\mathrm{\ on\ the\ hypersphere}")
    G = F("grad_norm")
    hm1 = heatmap!(ax1, αs, θs, G; colormap = :viridis, rasterize = 2)
    contour!(ax1, αs, θs, G; levels = 12, color = (:white, 0.5), linewidth = 0.6)
    mark_coalescences!(ax1; labels = false)
    Colorbar(fig[1, 2], hm1)
    ax2 = angle_axis(fig[1, 3]; title = L"\mathrm{(b)\ }\Lambda^2\psi_{3,0}\ \mathrm{(pseudo\!-\!log\ colour)}")
    hm2 = signed_map!(ax2, αs, θs, F("lambda2_psi30"); levels = nothing, colorscale = Makie.pseudolog10)
    mark_coalescences!(ax2; labels = false)
    Colorbar(fig[1, 4], hm2)
    save_tight("psi30_derivatives", fig)
end

function fig_recurrence()
    fig = Figure(size = (1150, 500))
    ax1 = angle_axis(fig[1, 1]; title = L"\mathrm{(a)\ }|(\Lambda^2-21)\psi_{3,0}-\mathrm{RHS}|")
    hm1 = log_map!(ax1, αs, θs, F("recurrence_residual_abs"); floor = 1e-16, colorrange = (1e-16, 1e-6))
    Colorbar(fig[1, 2], hm1; label = L"\mathrm{absolute\ residual}")
    ax2 = angle_axis(fig[1, 3]; title = L"\mathrm{(b)\ relative\ residual}")
    hm2 = log_map!(ax2, αs, θs, F("recurrence_residual_rel"); floor = 1e-17, colorrange = (1e-17, 1e-7))
    Colorbar(fig[1, 4], hm2)
    save_tight("psi30_recurrence", fig)
    r = F("recurrence_residual_abs")
    @printf("recurrence |residual|: median %.2e  max %.2e at %s\n", sort(vec(r))[end÷2],
        maximum(r), string((αs[argmax(r)[1]], θs[argmax(r)[2]])))
end

function fig_quadrature()
    fig = Figure(size = (620, 520))
    ax = angle_axis(fig[1, 1]; title = L"\mathrm{quadrature\ error\ estimate}")
    hm = log_map!(ax, αs, θs, F("psi30_quadrature_error"); floor = 1e-20, colormap = :viridis)
    Colorbar(fig[1, 2], hm)
    save_tight("psi30_quadrature_error", fig)
end

figures = isempty(ARGS) ? ["map", "hypersphere", "parts", "sectors", "sensitivity", "slices",
    "cusps", "zoom", "agreement", "derivatives", "recurrence", "quadrature"] : ARGS
for f in figures
    getfield(Main, Symbol("fig_", f))()
end
