# Admissibility and minimality after adding angular atoms

Codex, 2026-10-21. Supplements `codex/codex_minimality.md`.

**PROVED, conditional on reconstruction.** Let

    T A = exp(-sigma/2) sum_k log(sigma)^k sum_i g_i F_(k,i)

converge in a specified topology. For the same physical closed Coulomb
form q, define admissibility by

    T A in D(q),  q(phi,T A)=E <phi,T A> for all phi in D(q).

The representation theorem again identifies its image with ker(H-E).
At the simple ground eigenvalue that image is one dimensional. Adding
atoms changes the representation, not this statement. An equivalence to
“minimality” is valid only if that word denotes this admissibility
condition; no growing-block Pincherle theorem is gained.

**PROVED, gauge qualification.** The coefficient solution set need not be
one dimensional. Its kernel contains every allowed cancellation
sum_i g_i F_i=0, as well as the previous log-sector gauges. For example,
if the coefficient class admits g_L G, then F_one=-g_L G, F_L=G
reconstructs zero. Rational independence of a finite atom list only rules
out rational coefficient identities; it does not rule out identities
with unrestricted infinite Laguerre coefficient functions. One needs a
specified coefficient module, complement/normalization, and an injectivity
proof, or must work modulo ker(T). The local “log-free” radial condition
alone does not remove angular atom gauges.

**PROVED, interior versus boundary.** The rational connection and cleared
recurrence are identities on the open cone. The clearing factors vanish
on coordinate faces and at the corner. They cannot replace collision
conditions or exclude boundary-supported distributions. Individual U and
V and their analytic continuations have apparent singularities and
cancellations; the physical order-two combination is better behaved.
Require the *sum* to belong to the form domain, not each formal component
separately. Singlet exchange also acts on the atom functions and their
branches: imposing u<->v symmetry separately on every scalar coefficient
would require a compatible atom transformation matrix first.

**PROVED, conditional finite variational statement.** For any finite,
linearly independent trial set contained in D(q), the physical energy and
Gram matrices are Hermitian and G>0. Thus nonreal tail Schur complements
exist, including with growing block dimensions. If the trial functions
are dependent, first quotient/remove their Gram nullspace. No formal
triangular atom connection makes the physical Hermitian pencil triangular
or guarantees its shell bandwidth: angular-atom and log overlaps are
generally dense. Additional terms can improve approximation without
changing this distinction.

**CONDITIONAL / still open.** Strong resolvent convergence follows under
the same explicit nested form-core density hypothesis as before. The
new atoms may accelerate convergence of energies and Fock functionals;
they do not by themselves prove that hypothesis, global convergence of
the log tower, continuity of a21, or convergence of raw coefficient
Schur data. With six atoms (including the radial-log conversion atom)
the audit-system formal shell size becomes 6*d_N (5*d_N for the minimal
physical combination using Ltilde) before
exchange reduction and log sectors; rectangular growing couplings and
the infinite sector issue remain. Corner orders constrain reconstructed
jets, not a rule k<=N/2 for Laguerre shells.

**PROVED functional invariance.** The radial Mellin finite-part formula
for a21 depends on the reconstructed Psi and its projected source, so it
is invariant under atom gauges. For a trial state it remains an estimator
with a residual finite-part defect. Better energy alone still certifies
no digits of that corner quantity. This task creates no new atom-state
numerical a21 estimates; exported states are for the later cross-check.

## 2026-10-21 update after TASK_codex_4 / sync 4

The saved eight-atom recurrence now has independent actual-function checks:
48 actions under the distance Hamiltonian, maximum relative discrepancy
2.34e-39, and 144 exact connection-intertwining entries. The saved block
counts are 16/15/8, correcting the reported 17/12/8; its shell
half-bandwidth is 10. These results concern differential identities and
finite Laguerre actions on the open cone. They add no boundary condition
or asymptotic selection theorem.

The finite atoms2 Ritz sequence substantially improves approximation.
An independent driver using the same read-only weak-form assembly, finer
quadrature and higher precision reproduces N=11's error 2.43148e-17 and
positive physical Gram matrix. The independent Mellin functional on the
exported state gives 0.476748379513, about 5.0e-7 above the supplied
reference. This last result is numerical evidence for improved corner
information, without an a posteriori error bound. Both quantities refer
to finite trial states rather than a constructed minimal recurrence
solution. Their reproducibility does not prove form-core density,
coefficient convergence, or continuity of the finite-part functional in
the energy norm. A norm controlling its projected source and corner
remainders, or an explicit residual finite-part estimate, is still needed.

Observed log-error/log-size slopes near -8 support improved algebraic
behavior over the tested range; they establish no superalgebraic rate.
The Float64 psi30 variant's variational violations underscore the
hypothesis that the assembled values and gradients must represent the
same physical form-domain trial functions. Arbitrary symmetric numerical
matrices do not inherit the physical min–max bound.

No revision to the theoretical admissibility/gauge statements above is
required. The earlier sentence deferring atom-state estimates describes
TASK_codex_3; those estimates and their uncertainty have now been audited
in `REVIEW_4_codex.txt` and `sync4_a21.out`.

The further 512-bit h=1/32 check gives N=11 E-E0=2.4314823230015396e-17
(1.49e-27 from the h=1/24 result); its reconstructed state's finite-part
functional agrees through 15 printed decimals. This strengthens the finite
numerical evidence while leaving the theoretical qualifications unchanged.
