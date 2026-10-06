include("harmonic_proj.jl")
for ρ in (0.75, 0.9, 0.99), kind in (:inv, :log, :xlog)
    ω=1-ρ^2
    q=quadgk(γ->zonal(kind, ρ, ω, γ, 1-ρ*cos(γ))*sin(γ)*sin(3γ), 0, π; rtol=1e-13)[1]
    @printf("ρ=%.2f %-5s M closed % .15f  quad % .15f  diff %.1e\n", ρ, kind, Mzonal(kind, ρ, ω), q, Mzonal(kind, ρ, ω)-q)
end
