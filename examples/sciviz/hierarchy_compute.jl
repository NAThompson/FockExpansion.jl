# Computes and caches the Fock hierarchy for the hierarchy_*.jl figures:
#   output/hierarchy.nc          ψ_kp(α, θ) for k ≤ 12 from the spectral chain (4D stack,
#                                 raw and per-slice normalised), plus the chain's native
#                                 grid in the separating coordinates (X, Y) for k ≤ 24;
#   output/hierarchy_package.nc  the package's closed forms / quadratures for k ≤ 4: the
#                                 only copy of them in output/. The root holds them on a
#                                 91² grid (ψ40: 25²); /fock_grid on the 64² grid of
#                                 fock_truncated_*.jl, which read it from here.
# Run: julia -t 4 --project=examples/sciviz examples/sciviz/hierarchy_compute.jl
# Run this before fock_truncated_data.jl and fock_truncated_maps.jl.
include(joinpath(@__DIR__, "hierarchy_common.jl"))
include(joinpath(@__DIR__, "fock_truncated_lib.jl"))

const KMAX_STACK = 12
const KMAX_NATIVE = 24

function write_stack(path)
    g, ψ = he_chain(; n = 64, kmax = KMAX_NATIVE)
    αs = interior(181)
    θs = interior(181)
    ks = keys_upto(KMAX_STACK)
    t = @elapsed vals = interp_grid(g, ψ, ks, αs, θs)
    @printf("interpolated %d coefficients in %.1f s\n", length(ks), t)
    nk, np = KMAX_STACK + 1, KMAX_STACK ÷ 2 + 1
    stack = fill(NaN, length(αs), length(θs), np, nk)
    normed = fill(NaN, size(stack))
    for (k, p) in ks
        A = vals[(k, p)]
        stack[:, :, p+1, k+1] = A
        normed[:, :, p+1, k+1] = A ./ maximum(abs, A)
    end
    np_native = KMAX_NATIVE ÷ 2 + 1
    native = fill(NaN, g.n, g.n, np_native, KMAX_NATIVE + 1)
    for (k, p) in keys_upto(KMAX_NATIVE)
        native[:, :, p+1, k+1] = ψ[(k, p)]
    end
    isfile(path) && rm(path)
    NCDataset(
        path,
        "c";
        attrib = [
            "title" => "Fock hierarchy psi_kp(alpha, theta) for helium, spectral chain",
            "Z" => HE.Z,
            "E" => HE.E,
            "a21" => HE.a21,
            "free_harmonics" => "a21 set; harmonic parts of psi_k0 for even k >= 4 set to zero",
            "source" => "derivations/k4/fock_spectral.jl, n = 64 Chebyshev grid",
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
        defVar(
            ds,
            "p",
            collect(0:(np-1)),
            ("p",);
            attrib = ["long_name" => "power of log R"],
        )
        defVar(ds, "k", collect(0:(nk-1)), ("k",); attrib = ["long_name" => "power of R"])
        hint = view_hint(
            vis = "heatmap",
            scale = "symlog",
            colormap = "auto",
            symmetric = true,
        )
        # Static coalescence overlays on the (α, θ) stack and on the native (X, Y) grid.
        angular = write_coalescence!(ds)
        native_xy = write_coalescence!(ds, COALESCENCE_XY; columns = "X Y")
        defVar(
            ds,
            "psi",
            stack,
            ("alpha", "theta", "p", "k");
            attrib = field_attrib(
                stack,
                [
                    "long_name" => "ψ_{k,p}(α, θ) (NaN where p > k/2)",
                    hint...,
                    view_hint(
                        probe_dim = "k",
                        playback_dim = "k",
                        multiples_dim = "k p",
                    )...,
                ],
                angular;
                kind = :field,
            ),
        )
        defVar(
            ds,
            "psi_normalized",
            normed,
            ("alpha", "theta", "p", "k");
            attrib = field_attrib(
                normed,
                [
                    "long_name" => "ψ_{k,p}/max|ψ_{k,p}|",
                    view_hint(
                        vis = "heatmap",
                        scale = "linear",
                        colormap = "RdBu",
                        symmetric = true,
                        domain = [-1.0, 1.0],
                        probe_dim = "k",
                        playback_dim = "k",
                        multiples_dim = "k p",
                    )...,
                ],
                angular;
                kind = :field,
            ),
        )
        defVar(
            ds,
            "X",
            g.X,
            ("X",);
            attrib = ["units" => "rad", "long_name" => "X = α − β", PI_TICKS],
        )
        defVar(
            ds,
            "Y",
            g.Y,
            ("Y",);
            attrib = ["units" => "rad", "long_name" => "Y = α + β", PI_TICKS],
        )
        defVar(
            ds,
            "pn",
            collect(0:(np_native-1)),
            ("pn",);
            attrib = ["long_name" => "power of log R"],
        )
        defVar(
            ds,
            "kn",
            collect(0:KMAX_NATIVE),
            ("kn",);
            attrib = ["long_name" => "power of R"],
        )
        defVar(
            ds,
            "psi_native",
            native,
            ("X", "Y", "pn", "kn");
            attrib = field_attrib(
                native,
                [
                    "long_name" => "ψ_{k,p} on the chain's (X, Y) grid",
                    "comment" => "Chebyshev grid in separating coordinates X = α − β, Y = α + β (cos β = sin α cos θ); corners = coalescences",
                    hint...,
                    view_hint(
                        probe_dim = "kn",
                        playback_dim = "kn",
                        multiples_dim = "kn pn",
                    )...,
                ],
                native_xy;
                kind = :field,
            ),
        )
    end
    println("wrote ", path)
end

function write_package(path)
    pars = (Z = HE.Z, E = HE.E, a21 = HE.a21)
    αs = interior(91)
    θs = interior(91)
    fs = [
        "psi10" => (a, t) -> psi10(a, t; Z = pars.Z),
        "psi20" => (a, t) -> psi20(a, t; pars...),
        "psi21" => (a, t) -> psi21(a, t; Z = pars.Z),
        "psi30" => (a, t) -> psi30(a, t; pars..., rtol = 1e-13),
        "psi31" => (a, t) -> psi31(a, t; Z = pars.Z),
        "psi41" => (a, t) -> psi41(a, t; pars...),
        "psi42" => (a, t) -> psi42(a, t; Z = pars.Z),
    ]
    fields = Pair{String,Matrix{Float64}}[]
    for (name, f) in fs
        t = @elapsed A = grid(f, αs, θs)
        @printf("%s on %d points: %.1f s\n", name, length(A), t)
        push!(fields, name => A)
    end
    # ψ40 is ~0.1 s per point at step = 0.125; use a coarse grid, evaluated serially
    # (psi40 failed inside Threads.@threads).
    α40 = interior(25; ε = 0.05)
    # Points within ~0.02 of θ = π throw a DomainError (a quadrature node of the internal ψ30 call lands
    # on the boundary); those are stored as NaN and reported.
    function psi40_safe(a, b)
        try
            psi40(a, b; pars..., step = 0.125)
        catch err
            err isa DomainError || rethrow()
            @printf("psi40(%.4f, %.4f) failed: DomainError\n", a, b)
            NaN
        end
    end
    t = @elapsed P40 = [psi40_safe(a, b) for a in α40, b in α40]
    @printf("psi40 on %d points: %.1f s\n", length(P40), t)
    hint = Dict(
        name => [
            "long_name" => "ψ_{$(name[4]),$(name[5])}(α, θ), FockExpansion.jl";
            view_hint(vis = "heatmap", colormap = "auto")
        ] for (name, _) in fs
    )
    write_angular(
        path,
        αs,
        θs,
        fields;
        attribs = hint,
        global_attribs = [
            "title" => "FockExpansion.jl coefficients k <= 4 (helium ground state)",
            "Z" => HE.Z,
            "E" => HE.E,
            "a21" => HE.a21,
        ],
    )
    NCDataset(path, "a") do ds
        defVar(
            ds,
            "alpha40",
            α40,
            ("alpha40",);
            attrib = [
                "units" => "rad",
                "long_name" => "α (hyperangle), ψ_{4,0} grid",
                PI_TICKS,
            ],
        )
        defVar(
            ds,
            "theta40",
            α40,
            ("theta40",);
            attrib = [
                "units" => "rad",
                "long_name" => "θ (interelectron angle), ψ_{4,0} grid",
                PI_TICKS,
            ],
        )
        # The coalescence overlays again, in the coordinates of the coarse ψ40 grid.
        view40 = write_coalescence!(ds; columns = "alpha40 theta40", suffix = "_40")
        xyz40 = write_xyz!(ds, α40, α40, ("alpha40", "theta40"); suffix = "40")
        defVar(
            ds,
            "psi40",
            P40,
            ("alpha40", "theta40");
            attrib = field_attrib(
                P40,
                [
                    "long_name" => "ψ_{4,0} (a_{40} = a_{42} = 0, step = 0.125)",
                    "coordinates" => xyz40,
                    view_hint(vis = "heatmap", colormap = "auto")...,
                ],
                view40;
                kind = :field,
            ),
        )
    end
end

"""Append the group /fock_grid: ψ_{k,p}, k ≤ 4, on the (α, θ) grid of fock_truncated_*.jl,
evaluated exactly as `coefficients(α, θ; kmax = 4)` does."""
function write_fock_grid(path)
    αs, θs = fock_grid()
    C = Matrix{Vector{Vector{Float64}}}(undef, length(αs), length(θs))
    t = @elapsed Threads.@threads for j in eachindex(θs)
        for i in eachindex(αs)
            C[i, j] = coefficients(αs[i], θs[j]; kmax = 4)
        end
    end
    @printf("fock_grid coefficients on %d² points: %.1f s\n", length(αs), t)
    NCDataset(path, "a") do ds
        g = defGroup(
            ds,
            "fock_grid";
            attrib = [
                "comment" => "grid of fock_truncated_*.jl; theta stops at pi - 0.03 (psi40)",
            ],
        )
        # Own alpha/theta dimensions (64), shadowing the root ones (91).
        defDim(g, "alpha", length(αs))
        defDim(g, "theta", length(θs))
        defVar(
            g,
            "alpha",
            αs,
            ("alpha",);
            attrib = ["units" => "rad", "long_name" => "α (hyperangle)", PI_TICKS],
        )
        defVar(
            g,
            "theta",
            θs,
            ("theta",);
            attrib = ["units" => "rad", "long_name" => "θ (angle between r1, r2)", PI_TICKS],
        )
        xyz = write_xyz!(g, αs, θs, ("alpha", "theta"))
        angular = write_coalescence!(g; prefix = "/fock_grid")
        for (k, p) in FOCK_KP
            A = [c[k+1][p+1] for c in C]
            defVar(
                g,
                "psi$k$p",
                A,
                ("alpha", "theta");
                attrib = field_attrib(
                    A,
                    [
                        "long_name" => "ψ_{$k,$p}(α, θ), coefficient of R^{$k} log^{$p} R",
                        "coordinates" => xyz,
                        view_hint(vis = "heatmap", scale = "symlog", colormap = "auto")...,
                    ],
                    angular;
                    kind = :field,
                ),
            )
        end
    end
end

isfile(joinpath(OUT, "hierarchy.nc")) || write_stack(joinpath(OUT, "hierarchy.nc"))
write_package(joinpath(OUT, "hierarchy_package.nc"))
write_fock_grid(joinpath(OUT, "hierarchy_package.nc"))
