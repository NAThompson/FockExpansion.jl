# TASK 27 — exact flux route and accessible-orbit audit

2026-10-24. **Result: OPEN.** No minimal ODE of the actual `B_A`, and no
positive or negative H^1(E_q) identification, has been obtained. This is the
explicit-stall alternative in TASK 27. The obstruction is specified below;
it is not an ODE nonexistence result. All signs below use the positive
logarithm ell=log(-cos X), unless labelled `code`.

## 1. First variation versus the source proposed in the task

Put a=pi/2, b=3pi/2, I=(a,b), alpha=(X+Y)/2, beta=(Y-X)/2,
and M_z^(c)=d_z(cos z d_z)+(c/16)cos z. For the last electron-Coulomb
multiplication in ladder A,

    L4 phi4 := (M_X^(32)-M_Y^(32))phi4 = R4,
    R4 = (w v0/8) phi3 = (sqrt(2)/4)sin(alpha)cos(beta/2) phi3.

**PROVED:** this is the contribution of the specified multiplication to the
hierarchy; the complete physical `(4,0)` ladder equation can also contain
its `(4,1),(4,2)` hierarchy terms. These must either be retained in R4 or
proved zero/irrelevant in the selected path. `chain.jl` actually retains
them. In particular this report does not silently identify a multiplication
contribution with the entire resonant hierarchy.

Let Gk=N_b phik be the full locally analytic logarithmic variation, assuming
the degree-one property. Multiplication by wv0 commutes with N_b on the
chosen alpha,beta cover. Thus

    L4 G4 = (wv0/8) G3  [+ varied hierarchy terms if present],
    B_A(Y)=G4(b,Y),       g=G3=N_b phi3.

The task's source `-w*(-2v0*g)/16` is therefore a source for **G4**, not
for phi4. Moreover the saved A3 coefficient is only g(b,Y), not g(X,Y).
Even an exact ODE for that trace supplies neither the full interior g nor
the boundary value B_A of G4. This distinction is independent of the sign
convention: code coefficients are -B_A and -g(b,Y).

To see the local freedom directly, set t=X-b. Then cos X=t+O(t^3) and
-sin X=1+O(t^2). For an analytic G4=sum a_n(Y)t^n the coefficient recurrence
starts with

    a1 = M_Y^(32) a0 + R4^G(b,Y),     a0=B_A(Y).

At each subsequent order the new coefficient has nonzero multiplier
(n+1)^2. The equation determines normal Taylor coefficients *from* a0;
it does not determine a0 or its tangential ODE. The same issue occurs
for g, with c=21 and the varied order-two forcing. One must transport the
whole coupled varied hierarchy with its framing, or instead transport the
original phi3 and use the first-variation flux below.

**VERIFIED counterexample to substituting g for phi3:** take the analytic
test input phi3=sin(alpha), with no upper hierarchy terms. Its variation
g=0. For the sin(Y) mode of the order-four inverse analytic at a, however,
the positive-log coefficient is exactly

    16*(184-191*sqrt(2))/(1575*pi)
      = -0.27846287375427196017343286674345082327... != 0.

This is an exact test of the proposed source substitution, not a physical
ladder counterexample or a value of B_A. A separate 65-digit quadrature
checks the symbolic integral. A function's log variation alone loses the
analytic part that can generate a new log after multiplication/inversion.

## 2. Exact continuum flux and the missing kernel

Use the regular/log-free self-adjoint realization of M^(32) on I in
L^2(I,dY): the endpoint condition is vanishing limiting cos(Y) h'(Y).
It excludes the log solution. Write its real eigenfunctions h_j, eigenvalues
lambda_j and squared norms n_j=int_I h_j^2. Choose the same h_j(X) for
X in I. Each is regular at both endpoints, and h_j(b)!=0: the endpoint
regular power-series recurrence would otherwise force the zero solution.
The spectral expansion here is a continuum construction, not the finite
collocation eigensystem in `SOps`.

For any mode f_j of an actual continued solution,

    d_X{cos X*(h_j f_j' - h_j' f_j)} = h_j R_j.

If f_j is analytic at a, and has an integrable polylogarithmic source and
leading singular term B_j log(-cos X) at b, the endpoint fluxes are 0 and
h_j(b)B_j respectively. Terms t log^m t have vanishing flux. Hence

    B_j = int_a^b h_j(X) R_j(X) dX / h_j(b).

**PROVED under these explicit endpoint/source hypotheses.** For a general
start germ add its start flux divided by h_j(b). It cannot be omitted
when the input solution is already a varied branch with a log at a.
For reversed propagation b to a, cos'(a)=-1 contributes the corresponding
endpoint sign; one cannot reuse the orientation blindly.

Expanding R_j(X)=int_I h_j(eta)R(X,eta)deta/n_j gives the exact *formal*
two-dimensional source functional

    B(Y) = int_I int_I K32(Y;X,eta) R(X,eta) deta dX,
    K32 = sum_j h_j(Y)h_j(X)h_j(eta)/(n_j h_j(b)).

This is a precise spectral specification of the desired operator, subject
to convergence or distributional interpretation and interchange with the
source integral. It is **not** a closed evaluated Riemann kernel, a
convergence certificate, or an independent numerical algorithm. In
particular its appearance as a double integral does not solve task 1.
With the unvaried source it computes the first B; with the varied source
it computes a further log coefficient if the start-end hypotheses hold.
In the actual path-A saved chain the next variation is zero after
thresholding, rather than being B_A itself.

Three useful exact kernel modes at c=32 are

    h0=sin X, lambda0=0;
    h+=1+sqrt(2)cos X, lambda+=sqrt(2);
    h-=1-sqrt(2)cos X, lambda-=-sqrt(2).

They check normalization and give explicit projected fluxes. They do not
sum the infinite kernel. They also exhibit homogeneous ambiguity:
h_j(X)h_j(Y) solves L4 F=0 and is analytic at both X endpoints. Thus even
a known varied source together with analyticity at both X endpoints does
not determine G4. The physical monodromy lift supplies additional framing.
No inference here treats that ambiguity as permission to change the actual
physical solution.

**Precisely missing for the requested independent route:**

1. An evaluated Riemann/Poisson kernel for this regular endpoint-to-log
   map, including the corner conditions in Y, its sheet/path convention,
   and a convergence/error justification; K32 above specifies the map but
   does not evaluate it.
2. An independent exact or validated evaluator of the *unvaried* interior
   phi3^A (and any upper-log sources) on the selected strip. Alternatively,
   the full coupled transport of G3 and the framed lift to G4, which must
   determine its otherwise free a0=B_A. The boundary A3 ODE alone is not
   this transport or framing.
3. For either route, treatment of the alpha=pi,beta=0 crossing and a
   justified exchange of endpoint limits, spectral sums and source integrals.

No intrinsic impossibility is asserted: the original source history does
fix B_A. These are the unresolved constructions in this attempt.

## 3. Saved input and elliptic operator audit

Read-only inspection of `lad_B1.jls` confirms the four keys A,B,C,D; A
has 80 Chebyshev coefficients, stored at **256 bits** (about 77 decimal
digits), with last coefficient about -3.7563e-64. At ty=0,

    B_A^code = -0.00272518806189837190717990434226192011028214823379787...

This is a saved-series evaluation, not an independently certified value.
Increasing evaluation precision to 384 bits does not add source digits or
bound truncation error. The saved path-A chain has steps 0,1,2; its order-3
and order-4 values are zero at steps 1,2, after the program's 1e-50 cutoff.
`ladA_hi_160.log` is empty at inspection. No new high-precision coefficients
are inferred from it. Julia inspection peak RSS was 286,425,088 bytes.

`fac_A3.sage` rounds twice the floating coefficients of `op_A3.txt` and
then factors the resulting rational operator. Exact algebra on that rounded
operator is legitimate; without an identity or certified annihilation
argument it is not an exact certificate for the physical boundary g.
I therefore co-sign its Artin/Tate *model* interpretation conditionally,
not a new exact identification of the physical function. This does not
refute Claude's numerical factorization evidence.

For the x25 elliptic target

    E=4u^2(u^2-1)D^2+4u(3u^2-1)D+(4u^2+1),

**VERIFIED:** z=sqrt(u)*p transforms E p=0 into

    [4u^2(u^2-1)D^2+8u^3 D+u^2]z=0,

or, after division by u^2 on the generic arc,
`[4(u^2-1)D^2+8u D+1]z=0`. The choice of sqrt(u) branch on (-1,0)
is a constant complex normalization locally. An elliptic constituent can
also occur after a differential gauge or in a middle position, so direct
right-divisibility by the displayed scalar E is not a complete factor test.

**OPEN:** no minimal annihilator, no certified-gap guess from an independent
100-digit integral evaluator, and no differential-module factorization
of B_A. A bounded unsuccessful guess or the saved tail cannot replace
these. The report deliberately does not recycle Claude's ODE guesses as
an independent route.

## 4. What the accessible-corner orbit does and does not establish

Claude's statement concerns the physical solution and its four tested edge
meridians, together with finite deck images. Co-sign the observed small
orbit **NUMERICALLY**, subject to the sign and cutoff qualifications; do
not promote it to a theorem that these loops generate the complete
monodromy action on the cyclic differential module. Marching instability
on another sheet is an obstruction to that algorithm, not a topological
proof of inaccessibility by complex continuation.

For an exact version, let V_S be the span of all words in this specified
subgroup S applied to the physical germ (including every required deck
conjugate). If every generator's action on every displayed element is
certified, then the proposed finite span is invariant under S. Checking
only the alternating chain does not even establish that closure. Taking
one normal boundary trace loses further information, including positive
normal orders, as in x25's delta*period counterexample.

Conditional consequences are precise:

* If the entire cyclic continuation space under **all** monodromy were
  exactly three-dimensional, it could have no rank-five constituent in
  a submodule, middle, or quotient position. With deck images, the
  dimension must be counted; three symbolic names are not a dimension
  bound of three.
* If a complete, equivariant filtration under **all** monodromy had only
  finite/Tate factors (including its deck action), none of its subquotients
  could be T5. Unknown tangential structure of N_b Psi is part of what
  would need proving, not an available finite/Tate hypothesis.
* A terminating word on the selected physical germ gives a relation on
  that vector. To exclude a simple factor in any position, the tested
  relation must become a sufficiently restrictive operator identity on
  the entire cyclic space and remain restrictive on that factor. The
  selected vector can be in the word's kernel. x26's shared-line T5
  example already permits nonzero double and zero triple.
* For a T5 **submodule**, an actual differential extraction together with
  the applicable full physical analyticity/loop hypotheses can invoke
  the earlier submodule argument. This corner chain supplies no such
  extraction. A **quotient** or **middle** constituent requires the
  equivariant projection/subquotient and the physical pairing, or a
  full-space identity; absence from a tested scalar branch/trace is
  insufficient. D-module and solution-local-system arrows must be
  reversed consistently.

A simple reason not to confuse normal and tangential evidence is that
an elliptic family depending on u, tensored with a Tate log extension in
an independent normal variable, still has a short unipotent normal orbit.
Its elliptic tangential monodromy survives. This illustrates the logical
limitation; it is not a proposed physical decomposition.

**DISSENT** from an interpretation that the reported accessible orbit
excludes all T5 positions. **CO-SIGN** its value as numerical information
about those specific continued branches. **OPEN** whether B_A carries
H^1(E_q), and whether physical psi40 has a T5 constituent in any remaining
position. Conjecture A is not decided.

## Reproduction and scope

Serial commands, from the task root:

    PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 venv/bin/python codex/phase13/x27_exact.py
    PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 venv/bin/python codex/phase13/x27_quad.py
    OPENBLAS_NUM_THREADS=1 julia --startup-file=no --heap-size-hint=400M codex/phase13/x27_saved.jl
    PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase13/x27_audit.py

The first Julia inspection failed only on deserializing SGrid without its
definition; including the read-only spectral definitions fixed that cause.
No numerical chain was rerun and no input was edited. No agents, external
sources, package installation, commits or energy calculations were used.
Only x27_* outputs and append-only shared logs were written by this task.
