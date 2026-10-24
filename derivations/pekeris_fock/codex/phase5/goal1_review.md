# Goal 1 numerical review and handoff — 2026-10-22

PROVED result and precise limits: goal1_theorem.md. Literature/access: goal1_sources.md. Independent finite checks: goal1_verify.py and goal1_verification.md. This task does not pursue a21 digits or Goal 2's cyclic rank.

## VERIFIED independently (finite evidence)

The script reads only Codex's independently derived codex_recurrence.json (hash in the verification log), forms the Z=2, lambda=1, eps=1.704 physical polynomial pencil and checks cutoffs N=2,4,6 (7,22,50 singlet functions), serially at 45 decimal digits. It checks raw multiplicity-weighted Hermiticity, Gram positivity, H+5G positivity, tail positivity at the supplied E0, shell bandwidth 3, full backward Schur solves at c=-5 and E0+0.05i, normalization equivalence to deleting the first row, tail equation residuals, and decreasing ground constrained defect. Largest Schur equation residual <1.2e-45, direct-solve difference <7.8e-47. Peak RSS about 66 MiB. These are finite consistency checks, not numerical proof of the infinite theorem.

The raw dropped P row at N6 is -0.00018428754978695, independently reproducing Claude's -1.843e-4. At N6 and E0+0.05i it is -0.13982833921063+1.9004513494993i, reproducing Claude's -0.1398+1.900i. No Claude matrix/code is included in the independent script.

## VERIFIED by read-only review of Claude's contemporaneous results

Read claude5/goal1_backward.jl/.out and goal1_tail.jl/.out, not executed or modified. The backward file uses fixed eps, the general 57-term operator, exchange canonicalization, A000=1 and a deleted constant row, with shells above N set to zero. This is precisely the affine physical Galerkin problem in the theorem after multiplicity weighting. It is a square tail solve; it does not itself implement block-by-block elimination. The independent script verifies that proper Schur elimination gives the same finite answer.

At the snapshot backward.out contains N6..16 for all four energies. At E0 its raw residual magnitude decreases from 1.843e-4 to 5.217e-7. At E0+0.05i the final residual is approximately -0.1461+1.900i: finite low-shell convergence at a nonreal energy is fully compatible with a nonzero form-dual source. At real -2.5 it is about 68.28, with low-shell cutoff changes still substantial; no theorem of real off-spectrum tail convergence is inferred from that sample. At -1+0.5i the nonzero source likewise rules out an eigenfunction. The sample supports the criterion 'convergent physical normalized tail AND vanishing head residual', not 'small or decreasing shell coefficients alone'.

The tail snapshot through N18 compares perturbed exterior data, Miller(E0), Ritz vectors and roots of g_N. The printed g_N root and Ritz energy agree at every shown cutoff; at N18 both print -2.903724370803435. This is consistent with the exact Schur determinant identity near the ground branch. Shellwise Miller/Ritz discrepancies decrease in the shown sequence. Random exterior data are scaled by the shrinking last-shell norm, not by a fixed physical form norm; this is not an experiment establishing stability against all bounded tails or a uniform dominant/minimal ratio. The theorem instead proves stability when the exterior physical form norm tends to zero. No dominant-solution separation plot was present in this snapshot; none is fabricated.

## PROVED / OPEN task reconciliation

- PROVED: fixed-eps plain polynomial/Laguerre form core, using the supplied annihilator/Fourier argument rather than unread literature hypotheses.
- PROVED: each fixed raw Laguerre coefficient is form-continuous; the constant coefficient is nonzero on the positive ground state.
- PROVED: the normalized E0 tail form is coercive, hence the zero-tail finite-section/Miller criterion converges in physical form norm and fixed heads; its dropped residual vanishes. This is a substantive ground-state selection result, not only a renamed definition of admissibility.
- PROVED: nonreal physical resolvent/backward Schur convergence, ground Ritz/eigenprojection recovery, a best-approximation bound and a shifted inverse-iteration spectral rate with stated gap/overlap dependence.
- PROVED conditionally on actual finite atom admissibility, nesting and retention of the plain spaces: log/atom Galerkin form-core and physical resolvent convergence. This cannot be transferred automatically to the formal cleared sector recurrence.
- CITED with access limits: Klahn--Bingel I/II publisher abstracts, including Part II's explicit Theorem 10 locator; Klahn--Morgan DOI found but full text inaccessible. No unread exact exponent condition/rate used.
- OPEN: explicit shell-cutoff rate, arbitrary-formal-sequence Pincherle classification, ordinary limiting-series form convergence, and enriched recurrence gauge/physical-pencil identification. To close these respectively requires quantitative weighted approximation estimates, a reconstruction/dichotomy theorem, a weighted Laguerre basis theorem, and continuous injective enriched coordinates with admissible dense tails. Goal 2's order-three rank is outside TASK_codex_9.

All explicit writes were confined to codex/phase5/goal1_*.md/.py, NOTES_codex.md and append-only shared COLLAB.md / THEORY.md Log. Claude files and all external references were read only. No git metadata is available in this task directory; no git mutation, worktree, external write, or parallel process was used by Codex. Claude's own concurrent jobs updated its outputs independently.
