# Order-four reachability

2026-10-22; read-only source audit of src/feynman_k4.jl, fourth_order.jl,
and green_k4.jl. Distinguish finite holonomic representability from membership
in the particular order-three module.

**psi41^(2): yes in principle by the same finite direct-image method; no
membership in the existing order-three module is established.** The repository
now has a stronger input than a two-dimensional numerical Green integral:
feynman_k4.jl writes it as a one-dimensional integral of elementary functions.
The source has algebraic Feynman weight 1/sqrt(t(1-t)); rho, omega, c, ell
and the inverse trig arguments are algebraic after using coordinates
u=cos alpha and w=sin alpha cos theta. Cpow reduces finite trig powers;
Ktab/Jtab/zonal_I use only finite sums, log ell, angles and algebraic factors.
Differentiating the p-parameters finitely many times preserves this class.
Logarithms and inverse tangents of algebraic arguments are holonomic locally;
finite products, algebraic substitution and definite integrals preserve
holonomicity on nonsingular charts. Creative telescoping/relative de Rham
reduction consequently supplies a finite outer connection in principle.
The moving/singular endpoints and the resonant harmonic normalization must
be retained, rather than dropped in an integration-by-parts calculation.
No explicit reduced connection or minimality proof for this new integral has
been computed here. The repository's fitted PSI41_Z2_SHIFT is numerical
normalization data, not an exact derivation of its constants.

**psi40: similarly reachable in a larger finite holonomic extension in
principle, not solved by the 216-generator bound.** Its source is
12 psi41+2 psi42-2V psi30+2E psi20; its pure inverse is implemented using
the Green and iterated Green kernels. In the embedded S³ coordinates these
kernels are algebraic/trigonometric-angle expressions, with algebraic
coordinate measures. The azimuthal integrations and remaining angular
integrations are finite direct images of holonomic integrands once psi30
and the order-four source are supplied as holonomic functions. Singular
diagonal loci require a consistent integral/regularized continuation,
and finite-dimensional harmonic projections need endpoint/normalization
constants. This gives an existence route, not the explicit connection,
rank, clearing polynomial or stencil requested for an implemented atom system.
One cannot conclude “the same elliptic curve/alphabet suffices” from it:
new integration kernels can generate new periods/iterated integrals.

At resonance, a40 and a42 multiply the two singlet k=4 harmonics and are
independent global data. They are not determined by Z,E,a21 or the local
inhomogeneous recurrence. The source code explicitly takes a40,a42 as
arguments. Any claimed uniquely normalized psi40 must specify these two
numbers or use the zero-projection convention. psi41's harmonic projections
are fixed by next-order solvability, but the repository contains rounded
numerical pure-Z constants; their exact forms are not supplied by this step.

The same reasoning supports only fixed-order holonomicity. It does not prove
a common finite alphabet at all Fock orders, convergence of a global tower,
or the minimality/decay selection theorem discussed in earlier phases.
