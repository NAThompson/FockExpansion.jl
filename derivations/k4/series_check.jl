include("feynman41.jl")
for ρ in (0.1, 0.25, 0.29), kind in (:inv, :log, :xlog), γ in (0.4, 2.0)
    ω=1-ρ^2; c=cos(γ)
    cf=zonal(kind, ρ, ω, γ, 1-ρ*c)-(2/π)*Mzonal(kind, ρ, ω)*(4c^2-1)
    sr=zonal_series(kind, ρ, c)
    @printf("ρ=%.2f %-5s γ=%.1f closed % .15f series % .15f diff %.1e\n", ρ, kind, γ, cf, sr, cf-sr)
end
