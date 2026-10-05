# Exact harmonic (Y₄ₗ) content of the Feynman sum for ψ₄₁⁽²⁾.
# For a zonal U(p̂·x) = W(γ)/sin γ, Funk–Hecke gives ⟨Y, U⟩ = (4π/3) M Y(p̂),
# M = ∫₀^π W sin3γ dγ. By parts (I₁' = cos3γ r, I₂' = sin3γ r, r = r₀ + κ sin3γ, r₀ = -sinγ f/4):
#   M = (1/3)[(π/2)I₁(π) + κπ/8 - J₁/2 + J₂/12 + J₃/6],
#   J₁ = ∫γ cos3γ r₀ = π I₁⁰(π) - ∫₀^π I₁⁰ dγ,  J₂ = ∫sin6γ cos3γ r₀,  J₃ = ∫sin³3γ r₀,
# with I₁⁰ the κ-free I₁ = -(a1 log(1-ρ) + ac(c) log ℓ + poly(c))/4 (zonal_coeffs.py).
# Integrals over [0, π] of cos-polynomials against f use the Fourier coefficients of f.
include("feynman41.jl")

# (machinery moved into feynman41.jl)
