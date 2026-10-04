# FockExpansion.jl

This repo is my place to organize and validate new and old contributions to the Fock expansion for two electron atoms.

In 1954, Vladimir Fock proposed solving the Schrodinger equation for helium via an asymptotic expansion of the form

```math
\Psi(\rho, \alpha, \theta) = \sum_{k=0}^{\infty} \rho^k \sum_{p=0}^{\lfloor k/2 \rfloor} \ln(\rho)^{p} \psi_{k,p}(\alpha, \theta)
```

where $\rho := \sqrt{r_1^2+r_2^2}$, $\theta$ is the angle between electrons and $\alpha := 2 \mathrm{arctan}(r_2/r_1)$.
Later, Morgan demonstrated that this was not merely an asymptotic result, but was indeed pointwise convergent everywhere.
The "Fock coefficients" $\psi_{k,p}$ are related by a recurrence

```math
\left(\Lambda^2 - k(k+4)\right)\psi_{k,p} = 2(k+2)(p+1)\psi_{k,p+1} + (p+1)(p+2)\psi_{k,p+2} - 2V\psi_{k-1,p} + 2E\psi_{k-2,p}
```

where $\Lambda^2 = -4\Delta_{S^3}$ is the hyperspherical angular operator, $E$ is the energy, and $V = \rho\left(1/r_{12} - Z/r_1 - Z/r_2\right)$ depends only on the angles.
By taking $\psi_{0,0}(0,0) = 1$ to establish the global scale, we can then begin to recover analytic forms for all subsequent coefficients.

Fock was able to recover the $\psi_{1,0}$ and $\psi_{2,1}$ terms, at which point progress stalled.
In the 80s, Abbott and Maslen used an early version of Mathematica to recover $\psi_{2,0}$.
Recently, Liverts and coauthors managed to get $\psi_{3,1}$ as well as parts of $\psi_{3,0}$ computed.
They also give an integral formula found via Green's functions techniques to represent $\psi_{3,0}$ explicitly.
Langner's thesis gave $\psi_{3,0}$ as a doubly-infinite sum over terms involving ${}_3F_2$.

Recently I set `claude` and `codex` on the problem of finding a simple representation of $\psi_{3,0}$.
What we found cannot be considered _simple_, per se, but it is now represented as a finite sum of reasonably well-known special functions.
We also provide an numerical implementation in Julia which supports arbitrary precision and automatic differentiation.

Encouraged by these results, I figured I could probably press on beyond $\psi_{3,0}$ and perhaps slowly but surely turn the Schrodinger equation for helium into a special function.
We're far from that goal, but this is where we'll try to do it.


## $\psi_{3,0}$

The third-order coefficient is the sum of three pieces:

1. An explicit classical expression, using elementary functions and the Clausen function (a dilogarithm).
2. A short elementary expression containing the state parameters `E` and `a21`.
3. Four elliptic contributions.

The mathematical reduction expresses the elliptic contributions through one two-variable function and its derivative.
For numerical evaluation, this version uses the equivalent pair of moments `N₀` and `N₂`, combined into one integral per contribution.
The moments now use an equivalent elementary kernel with arctangent squares and finite endpoint limits, removing the logarithmic endpoint singularity.
Full evaluation takes about 14–15 microseconds.


## Available coefficients

| Function | Method |
|---|---|
| `psi00(α, θ)` | Constant 1 |
| `psi10(α, θ; Z)` | Elementary closed form |
| `psi20(α, θ; Z, E, a21)` | Classical closed form, using Clausen/dilogarithm identities |
| `psi21(α, θ; Z)` | Elementary closed form |
| `psi31(α, θ; Z)` | Elementary closed form |
| `psi30(α, θ; Z, E, a21)` | Explicit classical part plus four elliptic quadratures |
| `psi30_parts(α, θ; Z, E, a21)` | Separate contributions and estimated quadrature error |

The evaluators preserve Float64 or BigFloat input types and require `0 < α, θ < π`. Exact collinear and coalescence endpoints are not implemented. Near those boundaries, cancellation can reduce accuracy. `rtol` controls quadrature, not total relative error; the returned quadrature estimate is not a rigorous bound and excludes roundoff in the classical part. ForwardDiff gradients and Hessians are supported for `psi30`; quadrature controls the dual components as well as the value. For high precision, supply angles and state parameters as BigFloat and tighten `rtol`.

## Use

From this directory:

```sh
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

```julia
using FockExpansion

# Illustrative parameters, not a fitted ground state:
state = (Z=2.0, E=-2.9, a21=0.123)
psi20(0.8, 1.2; state...)
psi30(0.8, 1.2; state...)
psi30_parts(0.8, 1.2; state..., rtol=1e-12)
```

To use from another Julia environment:

```julia
using Pkg
Pkg.develop(path=expanduser("~/FockExpansion.jl"))
using FockExpansion
```

## Plots

Plotting dependencies live in a separate example environment:

```sh
julia --project=examples -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
julia --project=examples examples/angular_plots.jl
```

This creates PNG and SVG plots in `examples/output`. The example plots the lower coefficients and the three contributions to ψ₃₀ separately. Parameters are explicitly labelled as illustrative.

## Verification and provenance

The implementation translates the verified compact formulas in the research workspace's `round8_classical/compact_classical.py`, `round2/psi30_v2.py`, and the clean Python package's `general.py`. No Python or research-workspace dependency is needed at runtime.

Tests use three previously independently checked high-precision ψ₃₀ values, three 40-digit ψ₂₀ reference values, the exact Z=0 third-order formula, electron exchange symmetry, and tolerance refinement on a 25-point angular grid. Fourteen additional moment checks use independent 45-digit evaluations of the original logarithmic integral. Reference values are rounded only when stored as Float64. This is a translation audit, not a new proof of the underlying formula. The more ambitious function-class claims in the draft paper remain outside the package's guarantees.

Next additions: port the collinear classical-polylogarithm evaluator; add stable analytic endpoint limits; add further coefficients only with independent reference tests.

Verified locally with Julia 1.13.1: 128 assertions passed, and the CairoMakie example generated both PNG and SVG output. The PNG was visually inspected. Plot dependencies are pinned in the example Manifest; old cached plotting packages may be incompatible with Julia 1.13, so instantiate that environment before use.

### References

- V. A. Fock, On the Schrödinger equation of the helium atom, *Izv. Akad. Nauk SSSR Ser. Fiz.* **18**, 161 (1954).
- J. D. Morgan III, Convergence properties of Fock's expansion for S-state eigenfunctions of the helium atom, *Theor. Chim. Acta* **69**, 181 (1986).
- P. C. Abbott and E. N. Maslen, Coordinate systems and analytic expansions for three-body atomic wavefunctions: I. Partial summation for the Fock expansion in hyperspherical coordinates, *J. Phys. A* **20**, 2043 (1987).
- E. Z. Liverts and N. Barnea, Angular Fock coefficients: Refinement and further development, *Phys. Rev. A* **92**, 042512 (2015), [arXiv:1505.02351](https://arxiv.org/abs/1505.02351).
- E. Z. Liverts and N. Barnea, The Green's function approach to the Fock expansion calculations of two-electron atoms, *J. Phys. A* **51**, 085204 (2018), [arXiv:1705.09125](https://arxiv.org/abs/1705.09125).
- J. Langner, *Towards an exact solution to the Schrödinger Equation of the Helium atom*, Ph.D. thesis, National Yang Ming Chiao Tung University (2022).


## ForwardDiff recurrence residual

Run `julia --project examples/forwarddiff_residual.jl residual.csv` to apply
ForwardDiff's gradient and Hessian to the numerical `psi30` evaluator and compare
`(Λ²-21)ψ30` with `10ψ31-2Vψ20+2Eψ10`. The left side is differentiated directly;
it is not reconstructed from the recurrence. Scalar derivative rules for the
Clausen primitives are used, including removable limits.

Using the ground-state inputs `Z=2`, `E=-2.903724377034119598311159245194404446696925310`,
and `a21=0.47674787900`, the maximum absolute residual over three interior points
is 1.43e-14 (Float64), 9.41e-37 (128 bits), 4.08e-56 (192 bits), and 1.66e-75
(256 bits). This tests the angular recurrence for the supplied parameters;
it does not establish additional physical digits of the ground-state inputs.
The residual of a radially truncated Fock expansion is a separate quantity.

### Alternative third-order evaluators

`psi30_green(α, θ; Z, E, a21, rtol=1e-8)` evaluates the full Liverts
azimuthally integrated Green representation using Float64 nested quadrature.
It constructs the source `10ψ31 - 2Vψ20 + 2Eψ10` from the lower coefficients;
it does not call `psi30`. The integration uses exchange symmetry and splits
at the target point. `method=:duffy` is an alternative singularity transformation;
our initial timings favored the default `:rectangular` integration.
Requested quadrature tolerance is not the achieved error of the complete integral.

`table = LangnerTable(; Z=2.0)` loads the corrected double-series coefficients.
`psi30_langner(α, θ, table; E, a21, N=60, J=320)` evaluates a rectangular
truncation with vectorized cached powers. The table supports `N≤60, J≤320`;
convergence can be slow near boundaries, so increase truncations and compare
against a reference. `method=:horner` and `method=:direct` retain alternative
summation methods. `LangnerTable(; Z, T=BigFloat)` loads the stored decimal
coefficients at the caller's precision; the stored coefficients have finite
precision (100-digit construction, 90-digit storage).
The block convention and cubic sign in the printed thesis are corrected.

Run `examples/benchmark_representations.jl` for a warmed three-method comparison
against a high-precision reference. Float64 speed comparisons use Float64
accuracy targets; they do not establish 25-digit accuracy.
