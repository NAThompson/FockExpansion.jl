# Constant-field experiment (independent of claude6)

At (t,q)=(10,17),(13,22), over F_1004535809, all ten boundary
rows were reduced in the saved 45-row ordering. All eight independent
sign choices of sqrt(D),z1,z3 give identical 10x45 relations.

For pi=103,107,109,113,127, four samples determine rational fits and the
fifth is held out. At each point the 450 entries have the following
minimal regular training-fit degrees (numerator,denominator), and all
450 pass the held-out specialization:

* 10: (2,1)
* 40: (1,1)
* 400: (0,0)

The nonconstant denominators are exactly pi after monic normalization;
thus these fits have Laurent support contained in {-1,0,1}. These are
specialization tests, not a symbolic proof of pi dependence. No analytic
value of pi is inserted; pi is an independent parameter for reconstruction.

Conjugating i while keeping sqrt2 fixed leaves every entry invariant.
Conjugating sqrt2 negates precisely 50 entries and fixes the other 400;
conjugating both has the same effect. Hence there is no observed need
for i in the relation coefficients. The conservative field requested by
the task, Q(i,sqrt2)(pi)(t,q), is well posed. The experiments suggest the
smaller Q(sqrt2)(pi)(t,q); reconstruct the 50 odd entries after dividing
by the chosen sqrt2 embedding, and reconstruct pi Laurent components
separately. Do not silently treat pi=101 or 103 as an exact constant-field
reconstruction: neither detects the other Laurent components.

Changing the sqrt2 embedding at fixed t,q also changes the algebraic
coordinate parametrization of distances. This experiment concerns the
rational connection in the stated t,q chart, rather than holding physical
distances fixed while changing a numerical square-root convention.

All exact residues, fits, parity lists, basis/target indexing and test
records are in goal2x_independent.json and goal2x_constants.json. The
source order is the phase-five saved JSON, not a reordered intuitive
list. The represented rank/actual-function independence distinction
remains as in the phase-five handoff.
