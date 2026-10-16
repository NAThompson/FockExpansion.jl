# ILLUSTRATIVE patch at physical scale: the truncated Fock series F_3 blended with a
# global stand-in through a smooth radial cutoff, Z = χ F_3 + (1 - χ) W.
#
# W here is the Codex witness (outputs/a21_global/runs/o10_k8_d80_l1/result.json), the
# 151-term logarithmic Hylleraas trial that EFFECTIVE_SEED_CERTIFICATE certifies. It is
# NOT the heat-smoothed approximant of the construction (that cubature was never run),
# and the blend radius 0.1 bohr is chosen for visibility, not the certified r* = 2^-14.
#
#     julia -t 4 --project=examples/sciviz examples/sciviz/patched_blend.jl
include(joinpath(@__DIR__, "patched_common.jl"))

Wit = load_witness()
rb = 0.1                                          # χ = 1 on ρ ≤ rb/2, 0 on ρ ≥ rb
angles = [(π / 2, π / 2), (π / 4, π / 2), (π / 2, 2π / 3)]
anglabel = [
    L"\alpha=\pi/2,\ \theta=\pi/2",
    L"\alpha=\pi/4,\ \theta=\pi/2",
    L"\alpha=\pi/2,\ \theta=2\pi/3",
]
cols = [C.blue, C.orange, C.violet]

Rlog = 10 .^ range(-5, log10(1.5); length = 301)
Rlin = collect(range(1e-4, 1.5; length = 301))

"""Return W - 1, F_3 - 1 and Z - 1 along R at one angle (differences kept accurate)."""
function profile(Rs, α, θ; params = HE)
    c = [(k, p, f(α, θ)) for (k, p, f) in FOCK_TERMS]
    if params !== HE   # Fock coefficients with the witness's own E and a21
        c = [
            (1, 0, psi10(α, θ; Z = 2.0)),
            (2, 0, psi20(α, θ; params...)),
            (2, 1, psi21(α, θ; Z = 2.0)),
            (3, 0, psi30(α, θ; params..., rtol = 1e-12)),
            (3, 1, psi31(α, θ; Z = 2.0)),
        ]
    end
    w = zeros(length(Rs))
    Threads.@threads for i in eachindex(Rs)
        w[i] = Float64(witness_m1(Wit, Rs[i], α, θ))
    end
    f = [fock_m1(c, R, 3) for R in Rs]
    z = [chi(R, rb) * fi + (1 - chi(R, rb)) * wi for (R, fi, wi) in zip(Rs, f, w)]
    (; w, f, z)
end
lin = [profile(Rlin, a, t) for (a, t) in angles]
lg = [profile(Rlog, a, t) for (a, t) in angles]
wpar = (Z = 2.0, E = Float64(Wit.energy), a21 = Float64(Wit.a21))
lg_wpar = profile(Rlog, angles[1]...; params = wpar)
d = abs.(lg[1].f .- lg[1].w)
i = argmin(abs.(Rlog .- 2.0^-14))
@printf(
    "|F3 - W| at r* = 2^-14, angle 1: %.3e;  at 0.1: %.3e\n",
    d[i],
    d[argmin(abs.(Rlog .- 0.1))]
)

# ---------------------------------------------------------------- figure
fig = Figure(size = (1400, 950))
Label(
    fig[0, 1:3],
    "ILLUSTRATIVE — not a certified object: blend radius 0.1 bohr (certified r* = 2⁻¹⁴), global stand-in = Codex witness trial, not the heat-smoothed W";
    color = C.orange,
    fontsize = 15,
    font = :bold,
)
for (j, (a, t)) in enumerate(angles)
    ax = Axis(
        fig[1, j];
        xlabel = L"$R$ (bohr)",
        ylabel = j == 1 ? L"\Psi\ \ (\Psi(0)=1)" : "",
        title = anglabel[j],
    )
    vspan!(ax, rb / 2, rb; color = (C.yellow, 0.25))
    lines!(
        ax,
        Rlin,
        1 .+ lin[j].w;
        color = C.gray,
        linewidth = 5,
        label = L"global stand-in $W$ (witness)",
    )
    lines!(
        ax,
        Rlin,
        1 .+ lin[j].f;
        color = C.blue,
        linewidth = 2,
        linestyle = :dash,
        label = L"Fock $F_3$",
    )
    lines!(
        ax,
        Rlin,
        1 .+ lin[j].z;
        color = C.orange,
        linewidth = 2,
        label = L"blend $Z=\chi F_3+(1-\chi)W$",
    )
    ylims!(ax, 0, 1.05)
    j == 1 && axislegend(ax; position = :rt, framevisible = false, labelsize = 14)
    j == 1 && text!(
        ax,
        rb,
        0.03;
        text = "blend zone",
        align = (:left, :bottom),
        offset = (4, 0),
        color = C.ink2,
        fontsize = 13,
    )
end

ax = Axis(
    fig[2, 1:2];
    xscale = log10,
    yscale = log10,
    xlabel = L"$R$ (bohr)",
    ylabel = L"|F_3 - W|",
    title = L"mismatch in the overlap; $|Z - W| = \chi\,|F_3 - W|$ vanishes beyond $r_b$",
)
vspan!(ax, rb / 2, rb; color = (C.yellow, 0.25))
for (r, lab, y) in
    ((2.0^-14, L"r_* = 2^{-14}", 1e-1), (2.0^-13, L"2^{-13}\ \mathrm{(sharpened)}", 1e-17))
    vlines!(ax, [r]; color = C.gray, linestyle = :dash)
    text!(
        ax,
        r,
        y;
        text = lab,
        align = (:left, :top),
        offset = (4, 0),
        color = C.ink2,
        fontsize = 13,
    )
end
text!(
    ax,
    rb,
    1e-1;
    text = "blend zone",
    align = (:right, :top),
    offset = (-2, 0),
    color = C.ink2,
    fontsize = 13,
)
vspan!(ax, 1 / 4096, 1 / 2048; color = (C.aqua, 0.2))
text!(
    ax,
    1 / 4096,
    1e-1;
    text = "extraction\nshell",
    align = (:left, :top),
    offset = (2, -22),
    color = C.ink2,
    fontsize = 13,
)
for j in eachindex(angles)
    lines!(
        ax,
        Rlog,
        abs.(lg[j].f .- lg[j].w) .+ 1e-30;
        color = cols[j],
        linewidth = 2,
        label = anglabel[j],
    )
end
lines!(
    ax,
    Rlog,
    abs.(lg_wpar.f .- lg_wpar.w) .+ 1e-30;
    color = C.ink,
    linewidth = 2.5,
    linestyle = :dot,
    label = L"angle 1, $F_3$ with the witness's $E, a_{21}$ (coincides)",
)
text!(
    ax,
    1.5e-3,
    1e-13;
    text = L"$\approx 2\times10^{-5} R^2$ for $R \lesssim 10^{-3}$: the witness's finite",
    align = (:left, :bottom),
    color = C.ink2,
    fontsize = 13,
)
text!(
    ax,
    1.5e-3,
    1e-13;
    text = L"basis cannot reproduce $\psi_{2,0}$ exactly",
    align = (:left, :top),
    offset = (0, -2),
    color = C.ink2,
    fontsize = 13,
)
ylims!(ax, 1e-22, 3)
xlims!(ax, 1e-5, 1.5)
axislegend(ax; position = :rb, framevisible = false, labelsize = 13)

ax = Axis(
    fig[2, 3];
    xlabel = L"$R$ (bohr)",
    ylabel = "weight",
    title = L"radial cutoff $\chi$ ($C^\infty$)",
)
rr = range(0, 0.15; length = 400)
vspan!(ax, rb / 2, rb; color = (C.yellow, 0.25))
lines!(ax, rr, chi.(rr, rb); color = C.blue, linewidth = 2, label = L"\chi")
lines!(ax, rr, 1 .- chi.(rr, rb); color = C.gray, linewidth = 2, label = L"1-\chi")
axislegend(ax; position = :rc, framevisible = false)
rowsize!(fig.layout, 1, Relative(0.42))
savefig("patched_blend", fig)

# ---------------------------------------------------------------- 3D field for H5Web
αs = interior(33; ε = 1e-2)
θs = interior(33; ε = 1e-2)
R3 = 10 .^ range(-4, 0; length = 25)
ang = [(k, p, grid((a, t) -> f(a, t), αs, θs)) for (k, p, f) in FOCK_TERMS]
D = zeros(length(R3), length(αs), length(θs))
Zf = similar(D)
Threads.@threads for jt in eachindex(θs)
    for ia in eachindex(αs), ir in eachindex(R3)
        R = R3[ir]
        w = Float64(witness_m1(Wit, R, αs[ia], θs[jt]))
        f = sum(R^k * log(R)^p * A[ia, jt] for (k, p, A) in ang)
        D[ir, ia, jt] = log10(abs(f - w))  # -Inf where they agree exactly (NaN colour)
        Zf[ir, ia, jt] = 1 + chi(R, rb) * f + (1 - chi(R, rb)) * w
    end
end
ncgroup(
    "blend";
    attrib = [
        "comment" => "ILLUSTRATIVE only: F_3 (FockExpansion.jl, HE params) blended with the Codex witness trial at rb = 0.1 bohr; not the certified construction",
    ],
) do g
    putvar!(
        g,
        "R",
        Rlog,
        ("R",);
        attrib = ["units" => "bohr", "long_name" => "hyperradius R", axis_hint(Rlog)...],
    )
    putvar!(
        g,
        "angle",
        collect(1:3),
        ("angle",);
        attrib = [
            "long_name" => "angle pair (α, θ)",
            "comment" => "1:(pi/2,pi/2) 2:(pi/4,pi/2) 3:(pi/2,2pi/3)",
            "H5WEB_slice_labels" => "/blend/angle_label",
        ],
    )
    putvar!(
        g,
        "angle_label",
        ["(α, θ) = (π/2, π/2)", "(α, θ) = (π/4, π/2)", "(α, θ) = (π/2, 2π/3)"],
        ("angle",);
        attrib = ["long_name" => "angle pair label"],
    )
    # W and F3 are drawn with Z as series of one line plot.
    putvar!(
        g,
        "W",
        1 .+ hcat([p.w for p in lg]...)',
        ("angle", "R");
        attrib = [
            "long_name" => "Codex witness W, Ψ(0) = 1",
            "H5WEB_label" => "global stand-in W (witness)",
            "H5WEB_color" => "#8a8985",
        ],
    )
    putvar!(
        g,
        "F3",
        1 .+ hcat([p.f for p in lg]...)',
        ("angle", "R");
        attrib = [
            "long_name" => "Fock partial sum F_{3}, k ≤ 3",
            "H5WEB_label" => "Fock F_{3}",
            "H5WEB_color" => "#2a78d6",
            "H5WEB_dash" => Int8(1),
        ],
    )
    putvar!(
        g,
        "Z",
        1 .+ hcat([p.z for p in lg]...)',
        ("angle", "R");
        attrib = [
            "long_name" => "blend Z = χF_{3} + (1 − χ)W",
            "H5WEB_label" => "blend Z",
            "H5WEB_color" => "#eb6834",
            view_hint(vis = "line", x_dim = "R", series = "/blend/W /blend/F3")...,
        ],
    )
    # Reference marks on the mismatch plot (x = R): the blend zone, the extraction shell
    # and the Fock ball radius before and after sharpening.
    for (name, lohi, label, color) in (
        ("blend_zone", [rb / 2, rb], "blend zone", "#eda100"),
        ("extraction_shell", Float64.([A_SHELL / 2, A_SHELL]), "extraction shell", "#1baf7a"),
    )
        putvar!(
            g,
            name,
            lohi,
            ("end",);
            attrib = [
                "long_name" => "$label [lo, hi]",
                "units" => "bohr",
                "H5WEB_geometry" => "band",
                "H5WEB_axis" => "x",
                "H5WEB_label" => label,
                "H5WEB_color" => color,
            ],
        )
    end
    for (name, r, label) in (
        ("r_star", R_STAR, "r_{*} = 2^{-14}"),
        ("r_star_sharp", R_STAR_SHARP, "2^{-13} (sharpened)"),
    )
        putvar!(
            g,
            name,
            fill(Float64(r)),
            ();
            attrib = [
                "long_name" => "Fock ball radius $label",
                "units" => "bohr",
                "H5WEB_geometry" => "vline",
                "H5WEB_label" => label,
                "H5WEB_color" => "#8a8985",
                "H5WEB_dash" => Int8(1),
            ],
        )
    end
    putvar!(
        g,
        "abs_F3_minus_W",
        hcat([abs.(p.f .- p.w) for p in lg]...)',
        ("angle", "R");
        attrib = [
            "long_name" => "|F_{3} - W|",
            view_hint(
                vis = "line",
                scale = "log",
                x_dim = "R",
                multi_dim = "angle",
                overlays = "/blend/blend_zone /blend/extraction_shell /blend/r_star /blend/r_star_sharp",
            )...,
        ],
    )
    putvar!(
        g,
        "R3",
        R3,
        ("R3",);
        attrib = ["units" => "bohr", "long_name" => "hyperradius R", axis_hint(R3)...],
    )
    putvar!(
        g,
        "alpha",
        αs,
        ("alpha",);
        attrib = ["units" => "rad", "long_name" => "α (hyperangle)", PI_TICKS],
    )
    putvar!(
        g,
        "theta",
        θs,
        ("theta",);
        attrib = ["units" => "rad", "long_name" => "θ (interelectron angle)", PI_TICKS],
    )
    angular = write_coalescence!(g; prefix = "/blend")
    # Julia (alpha, theta, R3) is HDF5 (R3, theta, alpha): R3 is the slider.
    Dp = permutedims(D, (2, 3, 1))
    Zp = permutedims(Zf, (2, 3, 1))
    putvar!(
        g,
        "log10_abs_F3_minus_W",
        Dp,
        ("alpha", "theta", "R3");
        attrib = field_attrib(
            Dp,
            [
                "long_name" => "|F_{3} − W|, slider over R",
                view_hint(
                    vis = "heatmap",
                    colormap = "Inferno",
                    probe_dim = "R3",
                    playback_dim = "R3",
                    value_transform = "exp10",
                )...,
            ],
            angular;
            kind = :error,
        ),
    )
    putvar!(
        g,
        "Z_blend",
        Zp,
        ("alpha", "theta", "R3");
        attrib = field_attrib(
            Zp,
            [
                "long_name" => "blended Z (illustrative)",
                view_hint(probe_dim = "R3", playback_dim = "R3")...,
            ],
            angular;
            kind = :field,
        ),
    )
end
