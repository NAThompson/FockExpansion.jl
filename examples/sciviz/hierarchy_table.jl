# "Periodic table" of the Fock hierarchy: ψ_kp(α, θ) in a grid with rows k and columns p,
# each panel normalised by its own max |ψ_kp| (printed in the panel). Two versions:
# k ≤ 6 (hierarchy_table) and k ≤ 12 (hierarchy_table_k12).
# Needs output/hierarchy.nc from hierarchy_compute.jl.
# Run: julia -t 4 --project=examples/sciviz examples/sciviz/hierarchy_table.jl
include(joinpath(@__DIR__, "hierarchy_common.jl"))

function load_stack()
    NCDataset(joinpath(OUT, "hierarchy.nc")) do ds
        Array(ds["alpha"]), Array(ds["theta"]), Array(ds["psi"])
    end
end

"""LaTeX scientific notation with two significant digits."""
function sci(x)
    e = floor(Int, log10(abs(x)))
    -2 <= e <= 1 && return @sprintf("%.2g", x)
    @sprintf("%.1f\\times10^{%d}", x / 10.0^e, e)
end

function table_figure(αs, θs, stack, kmax; panel = 150)
    np = kmax ÷ 2 + 1
    fig = Figure(size = (panel * np + 160, panel * (kmax + 1) + 80), fontsize = 13)
    ticks = ([0, π / 2, π], [L"0", L"\pi/2", L"\pi"])
    for k = 0:kmax, p = 0:(k÷2)
        A = stack[:, :, p+1, k+1]
        m = maximum(abs, A)
        ax = Axis(
            fig[k+1, p+1];
            aspect = DataAspect(),
            xticks = ticks,
            yticks = ticks,
            xticklabelsvisible = k == kmax,
            yticklabelsvisible = p == 0,
            xticklabelsize = 11,
            yticklabelsize = 11,
            title = latexstring("\\psi_{$k,$p}\\;\\;", sci(m)),
            titlesize = 13,
            titlegap = 2,
            titlefont = :regular,
        )
        heatmap!(ax, αs, θs, A ./ m; colormap = DIVERGING, colorrange = (-1, 1), rasterize = 2)
        if minimum(A) < 0 < maximum(A)
            contour!(ax, αs, θs, A; levels = [0.0], color = :black, linewidth = 0.8)
        end
    end
    for k = 0:kmax
        Label(
            fig[k+1, 0],
            latexstring("k=$k");
            rotation = 0,
            tellheight = false,
            fontsize = 14,
        )
    end
    for p = 0:(np-1)
        Label(fig[0, p+1], latexstring("p=$p"); tellwidth = false, fontsize = 14)
    end
    # The upper-right triangle is empty: put the colour bar there.
    Colorbar(
        fig[1:3, np];
        colormap = DIVERGING,
        colorrange = (-1, 1),
        halign = :left,
        label = L"\psi_{k,p}/\max|\psi_{k,p}|",
        width = 14,
        tellheight = false,
        tellwidth = false,
    )
    Label(
        fig[kmax+2, 1:np],
        L"$\alpha$ horizontal, $\theta$ vertical. Panel title: $\psi_{k,p}$ and $\max|\psi_{k,p}|$. Black line: $\psi_{k,p}=0$.";
        fontsize = 13,
    )
    colgap!(fig.layout, 6)
    rowgap!(fig.layout, 6)
    fig
end

αs, θs, stack = load_stack()
savefig("hierarchy_table", table_figure(αs, θs, stack, 6; panel = 170))
savefig("hierarchy_table_k12", table_figure(αs, θs, stack, 12; panel = 110))
