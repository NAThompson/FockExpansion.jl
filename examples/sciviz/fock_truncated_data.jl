# Build the truncated series Ψ_K(R, α, θ), K = 0..4, on a 64³ grid for the H5Web viewer,
# from the angular coefficients ψ_{k,p} (k ≤ 4) on a 64 × 64 grid.
#
#     julia -t 4 --project=examples/sciviz examples/sciviz/fock_truncated_data.jl
#
# Needs output/hierarchy_package.nc from hierarchy_compute.jl (its /fock_grid group holds
# the ψ_{k,p}). Writes output/fock_truncated.nc.
include("fock_truncated_lib.jl")

const N = 64
(; αs, θs, coeffs) = load_fock_grid()
@assert (αs, θs) == fock_grid(N)

# Hyperradius: linear grid up to R = 1.5, where the series visibly degrades.
Rs = collect(range(1.5 / N, 1.5; length = N))
Ks = collect(0:4)

# Psi[iθ, iα, iR, iK] in Julia order appears in the file as (K, R, α, θ).
Psi = [psiK(coeffs[i, j], K, R) for j = 1:N, i = 1:N, R in Rs, K in Ks]
# Local energy error (H - E)Ψ_K / Ψ_K.
Eloc = [
    residual(coeffs[i, j], K, R, potential(αs[i], θs[j])) / psiK(coeffs[i, j], K, R) for
    j = 1:N, i = 1:N, R in Rs, K in Ks
]

path = joinpath(OUT, "fock_truncated.nc")
isfile(path) && rm(path)
NCDataset(
    path,
    "c";
    attrib = [
        "title" => "Truncated Fock expansion of the helium ground state",
        "Z" => Z,
        "E" => E,
        "a21" => A21,
        "a40" => 0.0,
        "a42" => 0.0,
        "theta_range" => "theta in [1e-3, pi - 0.03]: psi40 is not evaluable closer to pi",
        "coefficients" => "psi_kp on this grid: hierarchy_package.nc, group /fock_grid",
        "convention" =>
            "Psi_K = sum_{k<=K} R^k sum_p log(R)^p psi_{k,p}(alpha, theta); " *
            "r1 = R cos(alpha/2), r2 = R sin(alpha/2), psi00 = 1",
    ],
) do ds
    defVar(ds, "K", Ks, ("K",); attrib = ["long_name" => "truncation order K"])
    defVar(ds, "R", Rs, ("R",); attrib = ["units" => "bohr", "long_name" => "hyperradius R"])
    defVar(ds, "alpha", αs, ("alpha",); attrib = ["units" => "rad", "long_name" => "α (hyperangle)", PI_TICKS])
    defVar(ds, "theta", θs, ("theta",); attrib = ["units" => "rad", "long_name" => "θ (angle between r1, r2)", PI_TICKS])
    # The ψ_{k,p} themselves are in hierarchy_package.nc /fock_grid (one copy only).
    angular = write_coalescence!(ds)
    defVar(
        ds,
        "Psi",
        Psi,
        ("theta", "alpha", "R", "K");
        attrib = field_attrib(
            Psi,
            [
                "long_name" => "truncated Fock series Ψ_{K}(R, α, θ)",
                # One panel per truncation order K (heatmaps over α, θ), R on the slider.
                view_hint(
                    vis = "heatmap",
                    colormap = "Viridis",
                    probe_dim = "K",
                    playback_dim = "R",
                    multiples_dim = "K",
                    multiples_columns = 3,
                )...,
            ],
            angular;
            kind = :field,
        ),
    )
    defVar(
        ds,
        "local_energy_error",
        Eloc,
        ("theta", "alpha", "R", "K");
        attrib = field_attrib(
            Eloc,
            [
                "long_name" => "(H - E)Ψ_{K} / Ψ_{K}, exact from the recurrence",
                "units" => "hartree",
                view_hint(
                    vis = "heatmap",
                    scale = "symlog",
                    colormap = "auto",
                    probe_dim = "K",
                    playback_dim = "K",
                )...,
            ],
            angular;
            kind = :error,
        ),
    )
end
println("wrote ", path)
