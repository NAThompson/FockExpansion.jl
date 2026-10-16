# Agreement and convergence of the spectral chain:
#   (top) |chain − package| on a grid for every coefficient the package evaluates (k ≤ 4);
#   (bottom) spectral convergence of the Float64 chain in the grid size n, against a
#   128-bit chain at n = 64, at 40 random interior points.
# Needs output/hierarchy_package.nc from hierarchy_compute.jl.
# Run: julia -t 4 --project=examples/sciviz examples/sciviz/hierarchy_convergence.jl
include(joinpath(@__DIR__, "hierarchy_common.jl"))
using Random

const KCONV = 6

# --- chain vs package ---
pk = NCDataset(joinpath(OUT, "hierarchy_package.nc")) do ds
    (
        α = Array(ds["alpha"]),
        θ = Array(ds["theta"]),
        α40 = Array(ds["alpha40"]),
        θ40 = Array(ds["theta40"]),
        f = Dict(
            (k, p) => Array(ds["psi$k$p"]) for
            (k, p) in [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (4, 1), (4, 2), (4, 0)]
        ),
    )
end
g, ψ = he_chain(; n = 64, kmax = 4)
pkeys = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (4, 0), (4, 1), (4, 2)]
diffs = Dict{Tuple{Int,Int},Matrix{Float64}}()
for key in pkeys
    αs, θs = key == (4, 0) ? (pk.α40, pk.θ40) : (pk.α, pk.θ)
    c = interp_grid(g, ψ, [key], αs, θs)[key]
    diffs[key] = abs.(c .- pk.f[key])
    @printf(
        "ψ%d%d: max |chain - package| = %.1e (max |ψ| = %.2g)\n",
        key...,
        maximum(diffs[key]),
        maximum(abs, pk.f[key])
    )
end

# --- convergence in n ---
rng = MersenneTwister(1)
pts = [(0.05 + (π - 0.1) * rand(rng), 0.05 + (π - 0.1) * rand(rng)) for _ = 1:40]
ref = setprecision(BigFloat, 128) do
    gb, ψb = he_chain(BigFloat; n = 64, kmax = KCONV)
    x, bw = cheb_nodes(BigFloat, gb.n)
    Dict(
        key => [Float64(interp(x, bw, ψb[key], big(a), big(t))) for (a, t) in pts] for
        key in keys_upto(KCONV)
    )
end
ns = 8:4:64
err = Dict(key => Float64[] for key in keys_upto(KCONV))
for n in ns
    # At small n the resonant eigenvalue pairs are not resolved; the solver refuses.
    gn, ψn = try
        he_chain(; n, kmax = KCONV)
    catch e
        println("n = $n: ", sprint(showerror, e))
        foreach(key -> push!(err[key], NaN), keys_upto(KCONV))
        continue
    end
    x, bw = cheb_nodes(Float64, n)
    for key in keys_upto(KCONV)
        e = maximum(
            abs(interp(x, bw, ψn[key], a, t) - r) for ((a, t), r) in zip(pts, ref[key])
        )
        push!(err[key], e / maximum(abs, ref[key]))
    end
end
for key in keys_upto(KCONV)
    @printf(
        "ψ%d%d: rel. error n=16 %.1e  n=36 %.1e  n=64 %.1e\n",
        key...,
        err[key][3],
        err[key][7],
        err[key][end]
    )
end

# --- figure ---
fig = Figure(size = (1350, 900), fontsize = 15)
ticks = ([0, π / 2, π], [L"0", L"\pi/2", L"\pi"])
crange = (-17, -9)
top = fig[1, 1] = GridLayout()
for (i, key) in enumerate(pkeys)
    αs, θs = key == (4, 0) ? (pk.α40, pk.θ40) : (pk.α, pk.θ)
    ax = Axis(
        top[1, i];
        aspect = DataAspect(),
        xticks = ticks,
        yticks = ticks,
        yticklabelsvisible = i == 1,
        title = latexstring("\\psi_{$(key[1]),$(key[2])}"),
        xlabel = L"\alpha",
        ylabel = i == 1 ? L"\theta" : "",
        xticklabelsize = 12,
        yticklabelsize = 12,
    )
    heatmap!(
        ax,
        αs,
        θs,
        log10.(max.(diffs[key], 1e-18));
        colormap = :magma,
        colorrange = crange,
    )
end
Colorbar(
    top[1, length(pkeys)+1];
    colormap = :magma,
    colorrange = crange,
    label = L"\log_{10}|\psi^{\mathrm{chain}}-\psi^{\mathrm{pkg}}|",
)
Label(
    top[0, 1:length(pkeys)],
    L"\mathrm{(a)\ spectral\ chain\ }(n=64)\mathrm{\ vs\ FockExpansion.jl\ (}\psi_{4,0}\mathrm{\ on\ a\ }25\times25\mathrm{\ grid,\ step}=0.125)";
    fontsize = 15,
)
colgap!(top, 6)

bottom = fig[2, 1] = GridLayout()
ax1 = Axis(
    bottom[1, 1];
    yscale = log10,
    xlabel = L"\mathrm{grid\ size\ } n",
    ylabel = L"\mathrm{max\ relative\ error}",
    title = L"\mathrm{(b)}\ p=0",
)
ax2 = Axis(
    bottom[1, 2];
    yscale = log10,
    xlabel = L"\mathrm{grid\ size\ } n",
    title = L"\mathrm{(c)}\ p\geq1",
    yticklabelsvisible = false,
)
for k = 1:KCONV
    lines!(
        ax1,
        ns,
        max.(err[(k, 0)], 1e-17);
        color = SERIES[k],
        linewidth = 2,
        label = latexstring("\\psi_{$k,0}"),
    )
    scatter!(ax1, ns, max.(err[(k, 0)], 1e-17); color = SERIES[k], markersize = 7)
end
for (i, key) in enumerate([(2, 1), (3, 1), (4, 1), (4, 2), (5, 1), (6, 1), (6, 2), (6, 3)])
    lines!(
        ax2,
        ns,
        max.(err[key], 1e-17);
        color = SERIES[i],
        linewidth = 2,
        label = latexstring("\\psi_{$(key[1]),$(key[2])}"),
    )
end
linkyaxes!(ax1, ax2)
ylims!(ax1, 1e-16, 1)
axislegend(ax1; position = :rt, framevisible = false, nbanks = 2, labelsize = 13)
axislegend(ax2; position = :rt, framevisible = false, nbanks = 2, labelsize = 13)
Label(
    bottom[1, 3],
    "reference: 128-bit chain, n = 64; 40 random interior points";
    tellheight = false,
    fontsize = 13,
    rotation = π / 2,
)
rowsize!(fig.layout, 1, Relative(0.4))
savefig("hierarchy_convergence", fig)
