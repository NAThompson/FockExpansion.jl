# Green's-function representation of ψ₃₀ (FockExpansion.psi30_green): kernel, source,
# integrand, error maps and rtol convergence against the elliptic reference psi30.
#     julia -t 4 --project=examples/sciviz examples/sciviz/codex_green.jl
# Writes figures codex_green_*.{pdf,png} and CSV data to output/codex_green/.
include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "codex_util.jl"))
using Printf

const GDIR = joinpath(OUT, "codex_green")
mkpath(GDIR)
const FE = FockExpansion
const PAR = (Z = HE.Z, E = HE.E, a21 = HE.a21)

reference(α, θ) = psi30(α, θ; PAR..., rtol = 1e-13)
green(α, θ; rtol = 1e-8, method = :rectangular) = FE.psi30_green(α, θ; PAR..., rtol, method)

# Source without the measure sin²x sin t/(16π) that green_source folds in.
source(x, t) = FE.green_source(x, t, PAR...) * 16π / (sin(x)^2 * sin(t))
# What quadgk integrates for target (α, θ), on the folded panel x ∈ [0, π/2].
integrand(α, θ, x, t) =
    FE.green_source(x, t, PAR...) *
    (FE.green_kernel(α, θ, x, t) + FE.green_kernel(α, θ, π - x, t))

function fig_kernel()
    α, θ = 0.7, 0.9
    xs = range(0.005, π - 0.005; length = 401)
    K = [FE.green_kernel(α, θ, x, t) for x in xs, t in xs]
    fig = Figure(size = (1100, 450))
    ax = Axis(
        fig[1, 1];
        xlabel = L"source $\alpha'$",
        ylabel = L"source $\theta'$",
        title = L"kernel $G(\alpha,\theta;\alpha',\theta')$, target $(0.7, 0.9)$",
        aspect = DataAspect(),
    )
    lim = 12.0
    hm = heatmap!(ax, xs, xs, clamp.(K, -lim, lim); colormap = :balance, colorrange = (-lim, lim))
    contour!(ax, xs, xs, K; levels = [0.0], color = :black, linewidth = 1)
    contour!(ax, xs, xs, K; levels = 2:2:lim, color = (:white, 0.6), linewidth = 0.7)
    scatter!(ax, [α], [θ]; marker = :xcross, color = :black, markersize = 12)
    Colorbar(fig[1, 2], hm; label = L"G \;(\mathrm{clipped\ at\ } \pm 12)")
    # Radial cuts: G grows like -c log r, a straight line on a log-r axis.
    ax2 = Axis(
        fig[1, 3];
        xscale = log10,
        xlabel = L"distance $r$ from target",
        ylabel = L"G",
        title = L"\mathrm{logarithmic\ singularity}",
    )
    rs = 10.0 .^ range(-12, -0.5; length = 200)
    for (φ, c) in zip((0, π / 4, π / 2, 3π / 4), Makie.wong_colors())
        lines!(ax2, rs, [FE.green_kernel(α, θ, α + r * cos(φ), θ + r * sin(φ)) for r in rs]; color = c,
            label = latexstring("\\varphi = $(round(φ / π; digits = 2))\\pi"))
    end
    axislegend(ax2; position = :rt, framevisible = false)
    colsize!(fig.layout, 3, Relative(0.38))
    savefig("codex_green_kernel", fig)
    K
end

function fig_source()
    xs = range(0.004, π - 0.004; length = 401)
    S = [source(x, t) for x in xs, t in xs]
    fig = Figure(size = (650, 520))
    ax = Axis(
        fig[1, 1];
        xlabel = L"\alpha'",
        ylabel = L"\theta'",
        aspect = DataAspect(),
        title = L"source $10\psi_{31} - 2V\psi_{20} + 2E\psi_{10}$",
    )
    hm = heatmap!(ax, xs, xs, symlog.(S); colormap = :balance, colorrange = (-3, 3))
    contour!(ax, xs, xs, S; levels = [0.0], color = :black, linewidth = 1)
    contour!(ax, xs, xs, symlog.(S); levels = 1.3:0.1:2.5, color = (:white, 0.5), linewidth = 0.6)
    text!(ax, π / 2, 0.42; text = L"r_{12}=0", align = (:center, :bottom), color = :white, fontsize = 14)
    text!(ax, 0.06, π / 2; text = L"r_2=0", align = (:left, :center), rotation = π / 2, color = :white, fontsize = 14)
    text!(ax, π - 0.06, π / 2; text = L"r_1=0", align = (:right, :center), rotation = π / 2, color = :white, fontsize = 14)
    Colorbar(fig[1, 2], hm; label = L"\mathrm{sign}(S)\,\log_{10}(1+|S|)", ticks = symlog_ticks(3))
    savefig("codex_green_source", fig)
    S
end

const TARGETS = [("generic", 0.7, 0.9), ("near r12=0", 1.5, 0.15), ("near nucleus", 0.08, 2.0)]

function fig_integrand()
    xs = range(0.002, π / 2; length = 300)
    ts = range(0.002, π - 0.002; length = 500)
    fig = Figure(size = (1050, 500))
    hm = nothing
    for (k, (name, α, θ)) in enumerate(TARGETS)
        F = [integrand(α, θ, x, t) for x in xs, t in ts]
        ax = Axis(
            fig[1, k];
            xlabel = L"\alpha'",
            ylabel = k == 1 ? L"\theta'" : "",
            aspect = DataAspect(),
            title = latexstring("\\mathrm{$(replace(name, " " => "\\ "))}\\ (\\alpha,\\theta) = ($α, $θ)"),
        )
        hm = heatmap!(ax, xs, ts, log10.(abs.(F) .+ 1e-12); colormap = :viridis, colorrange = (-6, 1), rasterize = 2)
        contour!(ax, xs, ts, F; levels = [0.0], color = :white, linewidth = 0.8)
        hlines!(ax, [θ]; color = :red, linestyle = :dash, linewidth = 1)
        vlines!(ax, [α]; color = :red, linestyle = :dash, linewidth = 1)
        k > 1 && hideydecorations!(ax; grid = false)
    end
    Colorbar(fig[1, 4], hm; label = L"\log_{10}|\,S\,(G + G^{\ast})\,\sin^2\!\alpha' \sin\theta'/16\pi|")
    Label(fig[0, 1:3], "Integrand of psi30_green on the folded panel; dashed: the cuts of method = :rectangular, white: sign changes"; fontsize = 15)
    savefig("codex_green_integrand", fig)
end

timed(f) = (t = time_ns(); v = f(); (v, (time_ns() - t) / 1e9))

function error_maps(n = 41)
    αs = interior(n)
    θs = interior(n)
    ref = grid(reference, αs, θs)
    out = Dict{String,Matrix{Float64}}("reference" => ref)
    for m in (:rectangular, :duffy)
        V = similar(ref)
        T = similar(ref)
        Threads.@threads for j in eachindex(θs)
            for i in eachindex(αs)
                V[i, j], T[i, j] = timed(() -> green(αs[i], θs[j]; rtol = 1e-8, method = m))
            end
        end
        out["green_$m"] = V
        out["time_$m"] = T
    end
    write_matrix_csv(joinpath(GDIR, "maps.csv"), out)
    αs, θs, out
end

function fig_error_maps(αs, θs, M)
    fig = Figure(size = (1000, 860))
    for (k, m) in enumerate(("rectangular", "duffy"))
        err = abs.(M["green_$m"] .- M["reference"])
        ax = Axis(fig[1, k]; title = latexstring("|\\psi^{\\mathrm{Green}}_{30} - \\psi_{30}|,\\ \\mathrm{$m},\\ \\mathrm{rtol}=10^{-8}"),
            xlabel = L"\alpha", ylabel = k == 1 ? L"\theta" : "", aspect = DataAspect())
        hm = heatmap!(ax, αs, θs, log10.(err .+ 1e-17); colormap = :inferno, colorrange = (-14, -6))
        k == 2 && Colorbar(fig[1, 3], hm; label = L"\log_{10}\,\mathrm{error}")
        ax = Axis(fig[2, k]; title = latexstring("\\mathrm{time\\ per\\ point},\\ \\mathrm{$m}"),
            xlabel = L"\alpha", ylabel = k == 1 ? L"\theta" : "", aspect = DataAspect())
        hm = heatmap!(ax, αs, θs, log10.(M["time_$m"] .* 1e3); colormap = :cividis, colorrange = (0, 3))
        k == 2 && Colorbar(fig[2, 3], hm; label = L"\log_{10}\,\mathrm{ms}")
    end
    savefig("codex_green_error_maps", fig)
end

const PROBES = [("generic", 0.7, 0.9), ("equator", 1.5, 1.5), ("near_r12_0.10", π / 2 - 1e-9, 0.10),
    ("near_r12_0.03", π / 2 - 1e-9, 0.03), ("near_r1_0.05", 0.05, 1.2), ("collinear", 1.0, π - 0.02)]

function rtol_sweep()
    rtols = 10.0 .^ (-3:-1:-12)
    rows = Any[]
    for (name, α, θ) in PROBES, m in (:rectangular, :duffy)
        ref = reference(α, θ)
        tasks = map(rtols) do r
            Threads.@spawn begin
                green(α, θ; rtol = r, method = m)  # warm
                v, t = timed(() -> green(α, θ; rtol = r, method = m))
                (name, String(m), r, v, abs(v - ref), t)
            end
        end
        append!(rows, fetch.(tasks))
    end
    write_rows_csv(joinpath(GDIR, "rtol_sweep.csv"), ["point", "method", "rtol", "value", "error", "seconds"], rows)
    # Cost and accuracy of the elliptic evaluator itself, over the same probes.
    erows = Any[]
    for r in (1e-6, 1e-8, 1e-11), (name, α, θ) in PROBES
        psi30(α, θ; PAR..., rtol = r)
        t = @elapsed for _ = 1:1000
            psi30(α, θ; PAR..., rtol = r)
        end
        push!(erows, (name, r, abs(psi30(α, θ; PAR..., rtol = r) - reference(α, θ)), t / 1000))
    end
    write_rows_csv(joinpath(GDIR, "elliptic_timing.csv"), ["point", "rtol", "error", "seconds"], erows)
    rtols, rows
end

const RTOL_TICKS = (10.0 .^ (-3:-3:-12), [latexstring("10^{$k}") for k = -3:-3:-12])

function fig_rtol(rtols, rows)
    fig = Figure(size = (1100, 460))
    ax1 = Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = L"\mathrm{rtol}", ylabel = L"|\psi^{\mathrm{Green}}_{30} - \psi_{30}|",
        xreversed = true, xticks = RTOL_TICKS, title = L"\mathrm{accuracy}")
    ax2 = Axis(fig[1, 2]; xscale = log10, yscale = log10, xlabel = L"\mathrm{rtol}", ylabel = "seconds per point",
        xreversed = true, xticks = RTOL_TICKS, title = L"\mathrm{cost}")
    lines!(ax1, rtols, rtols; color = :gray, linestyle = :dot, label = "error = rtol")
    names = unique(first.(rows))
    for (k, name) in enumerate(names), (m, mk, ls) in (("rectangular", :circle, :solid), ("duffy", :utriangle, :dash))
        sel = filter(r -> r[1] == name && r[2] == m, rows)
        c = Makie.wong_colors()[k]
        e = max.([r[5] for r in sel], 1e-17)
        scatterlines!(ax1, [r[3] for r in sel], e; color = c, marker = mk, linestyle = ls,
            label = m == "rectangular" ? replace(name, "_" => " ") : nothing)
        scatterlines!(ax2, [r[3] for r in sel], [r[6] for r in sel]; color = c, marker = mk, linestyle = ls)
    end
    Legend(fig[1, 3], [ax1.scene.plots[1]; [LineElement(color = c) for c in Makie.wong_colors()[1:length(names)]];
            MarkerElement(marker = :circle, color = :black); MarkerElement(marker = :utriangle, color = :black)],
        ["error = rtol"; replace.(names, "_" => " "); ":rectangular"; ":duffy"]; framevisible = false)
    savefig("codex_green_rtol", fig)
end

fig_kernel()
fig_source()
fig_integrand()
αs, θs, M = error_maps()
fig_error_maps(αs, θs, M)
fig_rtol(rtol_sweep()...)
