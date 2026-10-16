# Radial behaviour of the truncated Fock series Ψ_K at fixed electron-triangle shapes:
# build-up of Ψ_K, order-by-order terms, increments (convergence), the local-energy
# residual (H - E)Ψ_K / Ψ_K, and the Kato cusp ratios.
#
#     julia -t 4 --project=examples/sciviz examples/sciviz/fock_truncated_lines.jl
#
# The residual is evaluated exactly from the recurrence (fock_truncated_lib.jl) and checked
# against an independent Cartesian Laplacian of Ψ_K by ForwardDiff (markers in the figure).
include("fock_truncated_lib.jl")
using ForwardDiff, LinearAlgebra

# Shapes of the electron triangle. Exact collinear and coalescence angles are not
# implemented, so the special shapes are approached to within 1e-3 rad (0.03 rad for θ → π,
# where psi40 is not evaluable closer).
const SHAPES = [
    (α = π / 2, θ = π / 2, label = L"r_1 = r_2,\ \theta = \pi/2"),
    (α = π / 2, θ = π - 0.03, label = L"near collinear: $\theta = \pi - 0.03$, $r_1 = r_2$"),
    (α = π / 2, θ = 1e-3, label = L"$r_{12} \to 0$ ($\theta \to 0$, $r_1 = r_2$)"),
    (α = 1e-3, θ = π / 2, label = L"$r_2 \to 0$ ($\alpha \to 0$)"),
]
const CS = [coefficients(s.α, s.θ; kmax = 4, step40 = 0.0625) for s in SHAPES]
const VS = [potential(s.α, s.θ) for s in SHAPES]
const Kmax = 4
Klabel(K) = L"K = %$K"

"""One horizontal legend below a 2 × 2 grid of panels, taken from the first panel."""
function shared_legend!(fig)
    ax = content(fig[1, 1])
    Legend(fig[3, 1:2], ax; orientation = :horizontal, framevisible = false, tellheight = true)
end

# ---------------------------------------------------------------------------------------
# 1. Ψ_K(R) at fixed shape, linear axes.
let
    fig = Figure(size = (1000, 640))
    Rs = range(0, 1.5; length = 400)[2:end]
    for (n, s) in enumerate(SHAPES)
        ax = Axis(
            fig[div(n - 1, 2)+1, mod(n - 1, 2)+1];
            title = s.label,
            xlabel = L"R",
            ylabel = L"\Psi_K(R)",
        )
        for K = 0:Kmax
            lines!(ax, Rs, [psiK(CS[n], K, R) for R in Rs]; color = KCOLORS[K+1], label = Klabel(K), linewidth = 2)
        end
        ylims!(ax, -0.6, 1.2)
    end
    shared_legend!(fig)
    savefig("fock_truncated_along_R", fig)
end

# ---------------------------------------------------------------------------------------
# 2. Order-by-order contributions |R^k (log R)^p ψ_{k,p}| on log-log axes.
let
    fig = Figure(size = (1000, 700))
    Rs = 10 .^ range(-4, log10(1.5); length = 300)
    for (n, s) in enumerate(SHAPES)
        ax = Axis(
            fig[div(n - 1, 2)+1, mod(n - 1, 2)+1];
            title = s.label,
            xlabel = L"R",
            ylabel = L"|R^k (\log R)^p \psi_{k,p}|",
            xscale = log10,
            yscale = log10,
        )
        for k = 1:Kmax, p = 0:(length(CS[n][k+1])-1)
            y = [abs(R^k * log(R)^p * CS[n][k+1][p+1]) for R in Rs]
            lines!(
                ax,
                Rs,
                max.(y, 1e-30);
                color = KCOLORS[k+1],
                linestyle = (:solid, :dash, :dot)[p+1],
                linewidth = 2,
                label = L"\psi_{%$k,%$p}",
            )
        end
        ylims!(ax, 1e-14, 10)
    end
    shared_legend!(fig)
    savefig("fock_truncated_terms", fig)
end

# ---------------------------------------------------------------------------------------
# 3. Increments |Ψ_K - Ψ_{K-1}| = |R^K Σ_p log^p R ψ_{K,p}| and the crossing radius where
# the newest increment stops being smaller than the previous one.
"""Smallest R in (Rlo, Rhi) where |term K| = |term K-1|, by bisection on a log grid."""
function crossing(c, K; Rlo = 1e-3, Rhi = 5.0)
    g(R) = log(abs(term(c, K, R))) - log(abs(term(c, K - 1, R)))
    Rs = 10 .^ range(log10(Rlo), log10(Rhi); length = 2000)
    for i = 1:(length(Rs)-1)
        if g(Rs[i]) < 0 && g(Rs[i+1]) ≥ 0
            a, b = Rs[i], Rs[i+1]
            for _ = 1:60
                m = sqrt(a * b)
                g(m) < 0 ? (a = m) : (b = m)
            end
            return sqrt(a * b)
        end
    end
    NaN
end
let
    fig = Figure(size = (1000, 700))
    Rs = 10 .^ range(-3, log10(3); length = 400)
    for (n, s) in enumerate(SHAPES)
        ax = Axis(
            fig[div(n - 1, 2)+1, mod(n - 1, 2)+1];
            title = s.label,
            xlabel = L"R",
            ylabel = L"|\Psi_K - \Psi_{K-1}|",
            xscale = log10,
            yscale = log10,
        )
        for K = 1:Kmax
            lines!(ax, Rs, [max(abs(term(CS[n], K, R)), 1e-30) for R in Rs]; color = KCOLORS[K+1], linewidth = 2, label = Klabel(K))
        end
        Rx = crossing(CS[n], 4)
        Rx3 = crossing(CS[n], 3)
        @printf("%-40s  crossing R(3|2) = %.3f   R(4|3) = %.3f\n", s.label, Rx3, Rx)
        isnan(Rx) || vlines!(ax, [Rx]; color = :gray40, linestyle = :dash)
        isnan(Rx) || text!(ax, Rx, 1e-11; text = L"R_\times = %$(round(Rx; digits = 2))", rotation = π / 2, align = (:left, :bottom), fontsize = 13, color = :gray30)
        ylims!(ax, 1e-12, 10)
    end
    shared_legend!(fig)
    savefig("fock_truncated_increments", fig)
end

# ---------------------------------------------------------------------------------------
# 4. Local energy residual (H - E)Ψ_K / Ψ_K, with an independent ForwardDiff check.
function cartesian_psi(x, K)
    a, b = x[1:3], x[4:6]
    r1, r2 = norm(a), norm(b)
    c = coefficients(2atan(r2, r1), acos(dot(a, b) / (r1 * r2)); kmax = K)
    psiK(c, K, sqrt(r1^2 + r2^2))
end
"""(H - E)Ψ_K / Ψ_K from the Cartesian Laplacian (6D ForwardDiff Hessian)."""
function cartesian_local_error(α, θ, R, K)
    r1, r2 = R * cos(α / 2), R * sin(α / 2)
    x = [r1, 0.0, 0.0, r2 * cos(θ), r2 * sin(θ), 0.0]
    ψ = cartesian_psi(x, K)
    lap = tr(ForwardDiff.hessian(y -> cartesian_psi(y, K), x))
    r12 = norm(x[1:3] - x[4:6])
    (-lap / 2 + (1 / r12 - Z / r1 - Z / r2 - E) * ψ) / ψ
end
let
    fig = Figure(size = (1000, 720))
    Rs = 10 .^ range(-4, log10(1.5); length = 300)
    Rcheck = [1e-3, 1e-2, 0.1, 0.5]
    worst = 0.0
    for (n, s) in enumerate(SHAPES)
        ax = Axis(
            fig[div(n - 1, 2)+1, mod(n - 1, 2)+1];
            title = s.label,
            xlabel = L"R",
            ylabel = L"|(H-E)\Psi_K / \Psi_K|",
            xscale = log10,
            yscale = log10,
        )
        for K = 0:Kmax
            y = [abs(residual(CS[n], K, R, VS[n]) / psiK(CS[n], K, R)) for R in Rs]
            lines!(ax, Rs, max.(y, 1e-30); color = KCOLORS[K+1], linewidth = 2, label = Klabel(K))
            if K ≤ 3 && s.α > 0.01 && s.θ > 0.01  # AD check away from the coalescences
                # At exactly α = π/2, psi20 switches to a limit formula whose ForwardDiff
                # Hessian misses the α-curvature, so the check is made 1e-6 away.
                αc = s.α + 1e-6
                cc = coefficients(αc, s.θ)
                yad = [abs(cartesian_local_error(αc, s.θ, R, K)) for R in Rcheck]
                yex = [abs(residual(cc, K, R, potential(αc, s.θ)) / psiK(cc, K, R)) for R in Rcheck]
                worst = max(worst, maximum(abs.(yad .- yex) ./ yex))
                scatter!(ax, Rcheck, yad; color = :black, marker = :circle, markersize = 9, strokewidth = 0)
            end
        end
        # Reference slopes R^{K-1} for K = 2, 3, 4 (up to logarithms).
        for K = 2:4
            R0 = 2e-4
            y0 = abs(residual(CS[n], K, R0, VS[n]) / psiK(CS[n], K, R0)) * 0.08
            lines!(ax, [R0, 20R0], [y0, y0 * 20.0^(K - 1)]; color = :gray50, linestyle = :dot)
            text!(ax, 20R0, y0 * 20.0^(K - 1); text = L"R^{%$(K-1)}", fontsize = 12, color = :gray40, align = (:left, :center))
        end
        ylims!(ax, 1e-12, 1e5)
    end
    shared_legend!(fig)
    @printf("max relative difference, recurrence residual vs Cartesian ForwardDiff: %.2e\n", worst)
    Label(fig[0, :], "markers: independent Cartesian Laplacian of Ψ_K (ForwardDiff), K ≤ 3"; fontsize = 13, color = :gray30)
    savefig("fock_truncated_residual", fig)
end

# ---------------------------------------------------------------------------------------
# 5. Kato cusp ratios. Electron-nucleus: average over directions of r1 of ∂Ψ/∂r1 at r1 = 0,
# divided by Ψ(r1 = 0), must equal -Z; electron-electron: same with r12 at fixed centre
# s = (r1 + r2)/2 must equal +1/2. For Ψ_K the error is O(R^K) (up to logarithms).
# Ψ averaged over the direction (Gauss-Legendre in the cosine) is a smooth function of the
# small distance u; value and slope at u = 0 come from a quartic fit at u = h, 2h, ..., 5h.
const GL_X = [-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640]
const GL_W = [0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891]

function fit_value_slope(us, ys)
    A = [u^j for u in us, j = 0:(length(us)-1)]
    c = A \ ys
    c[1], c[2]
end

"""Shape (R, α, θ) for electron 1 at distance u from the nucleus, direction cosine μ to r2."""
function en_shape(r2, u, μ)
    R = sqrt(u^2 + r2^2)
    (R, 2atan(r2, u), acos(μ))
end
"""Shape for the two electrons at s ± u/2, with direction cosine μ between u and s."""
function ee_shape(S, u, μ)
    r1 = sqrt(S^2 + u^2 / 4 + S * u * μ)
    r2 = sqrt(S^2 + u^2 / 4 - S * u * μ)
    cosθ = (S^2 - u^2 / 4) / (r1 * r2)
    (sqrt(r1^2 + r2^2), 2atan(r2, r1), acos(clamp(cosθ, -1, 1)))
end

function cusp_ratio(shape, d, Ks; hrel = 2e-3, step40 = 0.0625)
    h = hrel * d
    us = h .* (1:5)
    avg = zeros(length(Ks), length(us))
    for (j, u) in enumerate(us), (μ, w) in zip(GL_X, GL_W)
        R, α, θ = shape(d, u, μ)
        c = coefficients(α, θ; kmax = maximum(Ks), step40)
        for (i, K) in enumerate(Ks)
            avg[i, j] += w / 2 * psiK(c, K, R)
        end
    end
    [(v = fit_value_slope(us, avg[i, :]); v[2] / v[1]) for i in eachindex(Ks)]
end

let
    Ks = 1:4
    ds = 10 .^ range(-2.5, log10(0.8); length = 16)
    en = zeros(length(ds), length(Ks))
    ee = zeros(length(ds), length(Ks))
    Threads.@threads for i in eachindex(ds)
        en[i, :] = cusp_ratio(en_shape, ds[i], Ks)
        ee[i, :] = cusp_ratio(ee_shape, ds[i], Ks)
    end
    fig = Figure(size = (1000, 420))
    ax1 = Axis(
        fig[1, 1];
        title = L"electron–nucleus: $\langle\partial_{r_1}\Psi_K\rangle/\Psi_K \to -Z$ as $r_1 \to 0$",
        xlabel = L"r_2",
        ylabel = L"|\mathrm{ratio} + Z|",
        xscale = log10,
        yscale = log10,
    )
    ax2 = Axis(
        fig[1, 2];
        title = L"electron–electron: $\langle\partial_{r_{12}}\Psi_K\rangle/\Psi_K \to \frac{1}{2}$ as $r_{12} \to 0$",
        xlabel = L"|\mathbf{r}_1 + \mathbf{r}_2|/2",
        ylabel = L"|\mathrm{ratio} - \frac{1}{2}|",
        xscale = log10,
        yscale = log10,
    )
    for (i, K) in enumerate(Ks)
        scatterlines!(ax1, ds, max.(abs.(en[:, i] .+ Z), 1e-16); color = KCOLORS[K+1], label = Klabel(K), markersize = 7)
        scatterlines!(ax2, ds, max.(abs.(ee[:, i] .- 0.5), 1e-16); color = KCOLORS[K+1], label = Klabel(K), markersize = 7)
        @printf("cusp K=%d: e-n error at d=%.3g: %.3e, at d=%.3g: %.3e;  e-e: %.3e, %.3e\n", K, ds[1], abs(en[1, i] + Z), ds[end], abs(en[end, i] + Z), abs(ee[1, i] - 0.5), abs(ee[end, i] - 0.5))
    end
    axislegend(ax1; position = :rb, framevisible = false)
    savefig("fock_truncated_cusps", fig)
end
