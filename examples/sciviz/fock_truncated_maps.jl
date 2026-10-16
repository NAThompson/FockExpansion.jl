# Two-dimensional pictures of the truncated Fock series Ψ_K:
#   * Ψ_K and its local-energy error in the (r1, r2) plane at fixed θ (quarter disks),
#   * the local-energy error on hyperspherical shells R = const over (α, θ),
#   * Ψ_K painted on the hypersphere of electron-triangle shapes,
#   * the radius R_× where the fourth-order increment overtakes the third.
#
#     julia -t 4 --project=examples/sciviz examples/sciviz/fock_truncated_maps.jl
#
# Needs output/fock_truncated.nc from fock_truncated_data.jl and output/hierarchy_package.nc
# from hierarchy_compute.jl. Writes output/fock_truncated_maps.nc.
include("fock_truncated_lib.jl")

const DATA = joinpath(OUT, "fock_truncated.nc")
isfile(DATA) || error("run fock_truncated_data.jl first")

# ---------------------------------------------------------------------------------------
# Cached angular data. Arrays in the file are (theta, alpha[, R, K]) in Julia order.
ds = NCDataset(DATA)
αg = ds["alpha"][:]
θg = ds["theta"][:]
Rg = ds["R"][:]
Eloc = ds["local_energy_error"][:, :, :, :]
Psi = ds["Psi"][:, :, :, :]
close(ds)
# ψ_{k,p} on the same grid, from the package file; CG[iθ, iα].
pk = load_fock_grid()
@assert (pk.αs, pk.θs) == (αg, θg)
CG = permutedims(pk.coeffs)

# ---------------------------------------------------------------------------------------
# 1, 2. The (r1, r2) plane at fixed θ, evaluated on a polar (R, α) grid.
const ROWS = [
    (θ = 0.03, label = raw"\theta = 0.03"),
    (θ = π - 0.03, label = raw"\theta = \pi - 0.03"),
]
Rp = collect(range(0.0, 1.2; length = 121))[2:end]
αp = interior(129; ε = 2e-3)
polar = map(ROWS) do row
    c = Vector{Any}(undef, length(αp))
    Threads.@threads for i in eachindex(αp)
        c[i] = coefficients(αp[i], row.θ; kmax = 4)
    end
    V = potential.(αp, row.θ)
    X = [R * cos(a / 2) for R in Rp, a in αp]
    Y = [R * sin(a / 2) for R in Rp, a in αp]
    P = [psiK(c[i], K, R) for R in Rp, i in eachindex(αp), K = 0:4]
    Res = [residual(c[i], K, R, V[i]) / psiK(c[i], K, R) for R in Rp, i in eachindex(αp), K = 0:4]
    (; X, Y, P, Res)
end

function quarter_axis(pos, title)
    Axis(
        pos;
        title,
        xlabel = L"r_1",
        ylabel = L"r_2",
        aspect = DataAspect(),
        limits = (0, 1.2, 0, 1.2),
        xticks = [0, 0.5, 1],
        yticks = [0, 0.5, 1],
    )
end

let
    fig = Figure(size = (1150, 640))
    for (r, row) in enumerate(ROWS), K = 1:4
        ax = quarter_axis(fig[r, K], L"K = %$K,\ %$(row.label)")
        P = polar[r].P[:, :, K+1]
        surface!(ax, polar[r].X, polar[r].Y, P; colormap = :viridis, colorrange = (-0.5, 1), shading = NoShading, rasterize = 2)
        contour!(ax, polar[r].X, polar[r].Y, P; levels = -0.4:0.1:0.9, color = (:white, 0.6), linewidth = 0.8)
        contour!(ax, polar[r].X, polar[r].Y, P; levels = [0.0], color = :red, linewidth = 1.5)
        K > 1 && hideydecorations!(ax; grid = false)
        r == 1 && hidexdecorations!(ax; grid = false)
    end
    Colorbar(fig[1:2, 5]; colormap = :viridis, limits = (-0.5, 1), label = L"\Psi_K(r_1, r_2, \theta)")
    savefig("fock_truncated_r1r2", fig)
end

let
    fig = Figure(size = (1150, 640))
    for (r, row) in enumerate(ROWS), K = 1:4
        ax = quarter_axis(fig[r, K], L"K = %$K,\ %$(row.label)")
        L = log10.(abs.(polar[r].Res[:, :, K+1]))
        surface!(ax, polar[r].X, polar[r].Y, L; colormap = :magma, colorrange = (-6, 2), shading = NoShading, rasterize = 2)
        contour!(ax, polar[r].X, polar[r].Y, L; levels = -6:1:1, color = (:white, 0.5), linewidth = 0.8)
        K > 1 && hideydecorations!(ax; grid = false)
        r == 1 && hidexdecorations!(ax; grid = false)
    end
    Colorbar(fig[1:2, 5]; colormap = :magma, limits = (-6, 2), label = L"\log_{10} |(H-E)\Psi_K/\Psi_K|")
    savefig("fock_truncated_r1r2_residual", fig)
end

# ---------------------------------------------------------------------------------------
# 3. Local-energy error on hyperspherical shells, over (α, θ), from the cache.
let
    Rshow = [0.1, 0.5]
    fig = Figure(size = (1150, 600))
    lim = 2.0
    for (r, R0) in enumerate(Rshow), K = 1:4
        iR = argmin(abs.(Rg .- R0))
        ax = Axis(
            fig[r, K];
            title = L"K = %$K,\ R = %$(round(Rg[iR]; digits = 3))",
            xlabel = L"\alpha",
            ylabel = L"\theta",
            aspect = DataAspect(),
            xticks = ([0, π / 2, π], [L"0", L"\pi/2", L"\pi"]),
            yticks = ([0, π / 2, π], [L"0", L"\pi/2", L"\pi"]),
        )
        A = log10.(abs.(Eloc[:, :, iR, K+1]))  # (θ, α)
        heatmap!(ax, αg, θg, permutedims(A); colormap = :magma, colorrange = (-4, lim), interpolate = true)
        K > 1 && hideydecorations!(ax; ticks = false)
        r == 1 && hidexdecorations!(ax; ticks = false)
    end
    Colorbar(fig[1:2, 5]; colormap = :magma, limits = (-4, lim), label = L"\log_{10} |(H-E)\Psi_K/\Psi_K|")
    savefig("fock_truncated_shell_residual", fig)
end

# ---------------------------------------------------------------------------------------
# 4. Ψ_K painted on the hypersphere of shapes, x = r1/R, (y, z) = r2/R (cos θ, sin θ).
let
    X = [hypersphere_xyz(a, t)[1] for t in θg, a in αg]
    Y = [hypersphere_xyz(a, t)[2] for t in θg, a in αg]
    Zs = [hypersphere_xyz(a, t)[3] for t in θg, a in αg]
    fig = Figure(size = (1200, 680))
    function shell!(pos, vals, title, crange; labels = true)
        ax = Axis3(
            pos;
            title,
            aspect = :data,
            azimuth = -0.35π,
            elevation = 0.18π,
            xlabel = labels ? L"r_1/R" : "",
            ylabel = labels ? L"r_2\cos\theta/R" : "",
            zlabel = "",  # z = r2 sin θ / R; its label would spill into the next panel
            xticks = [0, 1],
            yticks = [-1, 0, 1],
            zticks = [0, 1],
        )
        surface!(ax, X, Y, Zs; color = vals, colormap = :viridis, colorrange = crange, shading = NoShading, rasterize = 2)
        # Coalescence points: r2 = 0 at (1, 0, 0) and r12 = 0 at (1, 1, 0)/√2.
        scatter!(ax, [Point3f(1, 0, 0), Point3f(1 / √2, 1 / √2, 0)]; color = [:red, :white], strokecolor = :black, strokewidth = 1, markersize = 10)
        ax
    end
    iR = argmin(abs.(Rg .- 0.5))
    for K = 1:4
        shell!(fig[1, K], Psi[:, :, iR, K+1], L"\Psi_%$K,\ R = %$(round(Rg[iR]; digits = 2))", (0.0, 0.8); labels = K == 1)
    end
    Colorbar(fig[1, 5]; colormap = :viridis, limits = (0.0, 0.8), label = L"\Psi_K")
    Rs4 = [0.25, 0.5, 0.75, 1.0]
    for (n, R0) in enumerate(Rs4)
        iR = argmin(abs.(Rg .- R0))
        shell!(fig[2, n], Psi[:, :, iR, 5] ./ maximum(Psi[:, :, iR, 5]), L"\Psi_4/\max \Psi_4,\ R = %$(round(Rg[iR]; digits = 2))", (0.3, 1.0); labels = n == 1)
    end
    Colorbar(fig[2, 5]; colormap = :viridis, limits = (0.3, 1.0), label = L"\Psi_4 / \max_{\mathrm{shell}} \Psi_4")
    savefig("fock_truncated_hypersphere", fig)
end

# ---------------------------------------------------------------------------------------
# 5. Radius R_× where |R^4 ψ_4| first exceeds |R^3 ψ_3| (fourth-order increment overtakes
# the third): an empirical, shape-resolved estimate of where the series stops improving.
"""Smallest R in (Rlo, Rhi) where |term K| ≥ |term K-1|."""
function crossing(c, K; Rlo = 1e-3, Rhi = 10.0, n = 600)
    g(R) = log(abs(term(c, K, R))) - log(abs(term(c, K - 1, R)))
    Rs = 10 .^ range(log10(Rlo), log10(Rhi); length = n)
    for i = 1:(length(Rs)-1)
        if g(Rs[i]) < 0 && g(Rs[i+1]) ≥ 0
            a, b = Rs[i], Rs[i+1]
            for _ = 1:50
                m = sqrt(a * b)
                g(m) < 0 ? (a = m) : (b = m)
            end
            return sqrt(a * b)
        end
    end
    NaN
end
Rx4 = [crossing(CG[j, i], 4) for j in eachindex(θg), i in eachindex(αg)]
Rx3 = [crossing(CG[j, i], 3) for j in eachindex(θg), i in eachindex(αg)]
for (name, A) in (("R_x(4|3)", Rx4), ("R_x(3|2)", Rx3))
    v = filter(isfinite, A)
    @printf("%s over the shape grid: min %.3f  median %.3f  max %.3f  (NaN: %d)\n", name, minimum(v), sort(v)[div(end, 2)], maximum(v), count(isnan, A))
end
let
    fig = Figure(size = (1000, 460))
    for (n, (A, t)) in enumerate(((Rx3, L"R_\times:\ |R^3\psi_3| = |R^2\psi_2|"), (Rx4, L"R_\times:\ |R^4\psi_4| = |R^3\psi_3|")))
        ax = Axis(
            fig[1, n];
            title = t,
            xlabel = L"\alpha",
            ylabel = L"\theta",
            aspect = DataAspect(),
            xticks = ([0, π / 2, π], [L"0", L"\pi/2", L"\pi"]),
            yticks = ([0, π / 2, π], [L"0", L"\pi/2", L"\pi"]),
        )
        heatmap!(ax, αg, θg, permutedims(A); colormap = :cividis, colorrange = (0.5, 2.5), interpolate = true)
        contour!(ax, αg, θg, permutedims(A); levels = 0.75:0.25:2.25, color = :white, linewidth = 0.8, labels = true, labelsize = 11, labelcolor = :white)
    end
    Colorbar(fig[1, 3]; colormap = :cividis, limits = (0.5, 2.5), label = L"R_\times")
    savefig("fock_truncated_crossing_radius", fig)
end

# ---------------------------------------------------------------------------------------
# NetCDF for the viewer: the (r1, r2) quarter disks as parametric surfaces, and R_×.
path = joinpath(OUT, "fock_truncated_maps.nc")
isfile(path) && rm(path)
NCDataset(path, "c"; attrib = ["title" => "Truncated Fock series: (r1, r2) planes and crossing radius"]) do nc
    defVar(nc, "K", collect(0:4), ("K",); attrib = ["long_name" => "truncation order K"])
    defVar(nc, "R", Rp, ("R",); attrib = ["units" => "bohr", "long_name" => "hyperradius R"])
    defVar(nc, "alpha", αp, ("alpha",); attrib = ["units" => "rad", "long_name" => "α (hyperangle)", PI_TICKS])
    defVar(nc, "alpha_grid", αg, ("alpha_grid",); attrib = ["units" => "rad", "long_name" => "α", PI_TICKS])
    defVar(nc, "theta_grid", θg, ("theta_grid",); attrib = ["units" => "rad", "long_name" => "θ", PI_TICKS])
    angular = write_coalescence!(nc; columns = "alpha_grid theta_grid", suffix = "_grid")
    for (r, row) in enumerate(ROWS)
        tag = r == 1 ? "theta_small" : "theta_near_pi"
        θtext = r == 1 ? "0.03" : "π − 0.03"
        defVar(nc, "r1_$tag", polar[r].X, ("R", "alpha"); attrib = ["units" => "bohr", "long_name" => "r_{1}"])
        defVar(nc, "r2_$tag", polar[r].Y, ("R", "alpha"); attrib = ["units" => "bohr", "long_name" => "r_{2}"])
        P = permutedims(polar[r].P, (2, 1, 3))
        defVar(
            nc,
            "Psi_$tag",
            P,
            ("alpha", "R", "K");
            attrib = field_attrib(
                P,
                [
                    "long_name" => "Ψ_{K} at θ = $θtext, polar grid (R, α)",
                    "theta" => row.θ,
                    view_hint(
                        vis = "heatmap",
                        colormap = "Viridis",
                        probe_dim = "K",
                        playback_dim = "K",
                    )...,
                ];
                kind = :field,
            ),
        )
        Res = permutedims(polar[r].Res, (2, 1, 3))
        defVar(
            nc,
            "local_energy_error_$tag",
            Res,
            ("alpha", "R", "K");
            attrib = field_attrib(
                Res,
                [
                    "long_name" => "(H - E)Ψ_{K}/Ψ_{K} at θ = $θtext",
                    "theta" => row.θ,
                    view_hint(
                        vis = "heatmap",
                        scale = "symlog",
                        colormap = "auto",
                        probe_dim = "K",
                        playback_dim = "K",
                    )...,
                ];
                kind = :error,
            ),
        )
    end
    for (name, A) in (("crossing_radius_4_3", Rx4), ("crossing_radius_3_2", Rx3))
        defVar(
            nc,
            name,
            A,
            ("theta_grid", "alpha_grid");
            attrib = field_attrib(
                A,
                [
                    "long_name" => "R where order $(name[end-2]) passes order $(name[end])",
                    "comment" => "smallest R where the order-K increment exceeds the order-(K-1) one",
                    "units" => "bohr",
                    view_hint(vis = "heatmap", colormap = "Cividis")...,
                ],
                angular;
                kind = :field,
            ),
        )
    end
end
println("wrote ", path)
