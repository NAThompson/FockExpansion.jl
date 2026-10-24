# Brevity memo: a five-page letter

Use THEORY.md's phase-8 results and latest qualifications, not its superseded
“PARTIAL” paragraph or ANSATZ.md's earlier open density status.
The paper has two linked achievements: an explicit structured finite recurrence,
and a physical weak-pencil selection theorem. Their full identification is OPEN.
Suggested title: “Finite recurrences with Fock structure and form-domain selection
for helium.” Avoid promising that the formal order-3 recurrence already selects
exactly the physical state.

## 1. Drop entirely from the letter

- The cyclic-vector reconstruction story: rank-45 experiments, degree >350,
  radical-field bounds, failed interpolations, and abandoned sector allocations.
- A separate “No exponential dichotomy” section and tail-propagation plot.
  Keep one sentence: Fock nonanalyticity excludes geometric decay of the full
  plain Laguerre array; no basis-independent dichotomy theorem is established.
- The a21 estimate: it is a second observable with conventions and error questions
  unrelated to the main construction/selection result.
- Claims of a new convergence exponent, “physical iff algebraic decay,” or a
  major order-3 speedup. THEORY.md reports marginal per-size gains, no slope change.
- Repeated contribution bullets, wish lists, and an inventory of all open problems.
  End with only the formal-recurrence identification and quantitative rate gaps.
- AI workflow chronology, co-sign histories, timings, and cache sizes. Acknowledge
  assistance once; preserve mathematical attribution in the research record.

## 2. Move to supplement or repository

**Supplement:** precise log-free sector/gauge definition and convergence assumptions;
conditional finite-log-tower argument; order-2 minimal-rank proof and atom formulas;
fixed-endpoint involution proof; admissibility estimates for r12 V and delta^4;
form-core/Hardy/coercivity proofs; head normalization and physical-dependency quotient;
nonreal Schur/resolvent extension; corrected Ritz comparisons versus basis dimension.
Include enough proof detail here that the main theorem is assessable.

**Repository:** full rational connections, row clearings, stencils/generator,
coefficient comparisons, actual-atom checks, per-block counts, CRT/finite-field
experiments, and exploratory plots. Point to claude8/stencil3.jl and the exact
exports/checks in codex/phase8/. Prior boundary matching is outside this paper.

**Main text:** one sector equation; the structured construction table below;
three short propositions with explicit scope; one paragraph explaining the selector.
Do not reproduce the elementary 33/177/117/57 offset table as a second table.

## 3. Compact propositions and honesty labels

**Proposition 1 — finite structured construction [PROVED / VERIFIED].**
The order-2 differential module has minimal rank five. The order-3 classical
and four fixed-endpoint panel families admit explicit global rational connections
on the strict triangle; each panel has a sufficient 23-reachable-state presentation
[PROVED; no minimal-rank claim]. Their cleared log-sector Laguerre operators have
2,577,274 offsets and shell half-bandwidth at most 16 [VERIFIED]. Actual-atom
operator checks cover all five blocks; independent symbolic certification covers
selected blocks, not every panel stencil. This is finite spatial coupling at each
log level, not a finite-dimensional exact eigenproblem or a finite exact log tower.

**Proposition 2 — admissible structured enrichment [PROVED; OPEN limitation].**
Corrected order-2 products using r12 V and the classical order-3 atoms are
Coulomb-form admissible. Multiplication by delta^4 suffices for all four rational
panel families. Universal powers below one fail; the least sufficient integer
power is OPEN in {1,2,3,4}. Nested symmetry-correct physical spaces retaining the
plain polynomial core are form-dense after quotienting physical dependencies.
Exact unweighted psi30 angular-jet inclusion in a finite polynomial span of the
weighted atoms remains OPEN; form closure does not establish it. Weighting requires
the d(log w) connection correction, so do not silently reuse unweighted stencils.

**Theorem 3 — ground-state zero-tail selection [PROVED under cited spectral input].**
Fix eps>0 and an isolated simple positive ground state psi0 of energy E0. Plain
exponential-polynomial spaces form a Coulomb form core; ell=A000 is form-bounded
and ell(psi0) is nonzero. The physical weak-pencil solution at E0 with ell=1,
zero tail above shell N, and all nonhead weak rows converges in form norm to
psi0/ell(psi0), with a Cea best-approximation bound. The omitted-row residual
vanishes; the nearby ground-state Ritz zeros converge monotonically to E0.
The result extends to nested admissible form-dense enrichments with a bounded
head nonzero on psi0 (use an overlap head, not an assumed enriched A000).
Identification with the formal row-cleared enriched recurrence, a classification
of all formal solutions, and an explicit shell-cutoff rate remain OPEN.

No CONJ claim is needed in the letter. If mentioned, all-order nonvanishing of
Fock log powers is CONJ, not a consequence of the allowed Fock hierarchy.

## 4. The single visual: one construction table

Use this table beside Proposition 1; omit convergence and tail-separation figures.
It makes the finite order-3 construction concrete without selling an unproved rate.
Counts refer to the recorded unweighted operators, not corrected weighted stencils.

| Block | Sectors | Offsets | Shell half-bandwidth |
|---|---:|---:|---:|
| Classical | 21 | 215,406 | <=12 |
| Panels 1,2 (each) | 23 | 562,570 | <=16 |
| Panels 3,4 (each) | 23 | 618,364 | <=16 |
| Total | 113 allocated | 2,577,274 | <=16 |

Caption: sufficient presentations, not independent minimal ranks; row factors have
maximum degree 12. All-block actual-atom checks give relative discrepancies about
10^-30 to 10^-54; selected blocks have independent exact certification.
The total counts these five blocks, not a newly audited combined order-2/3 system.

## 5. Page budget (five pages including references)

| Section/content | Pages |
|---|---:|
| Title, abstract (90–110 words), introduction | 0.65 |
| Ansatz, one sector equation, conventions | 0.55 |
| Structured construction, table, domain qualification | 1.25 |
| Selection theorem and short coercivity/density proof sketch | 1.35 |
| Verification paragraph, limitations, outlook | 0.45 |
| Acknowledgements and references | 0.75 |
| Total | 5.00 |

For four pages, shorten construction/proof sketches and references; for six,
expand the physical-pencil versus formal-recurrence distinction and proof sketch.
Do not restore the reconstruction chronology or a separate numerics section.

## 6. Corrections to OUTLINE.md before drafting

- Title/abstract conflate the explicit recurrence with the proved physical selector.
  State the two results separately; their formal enriched identification is OPEN.
- “Built-in decay” overstates the exponential gauge: unrestricted sums can undo it.
  Physical form-domain admissibility supplies the boundary condition.
- Section 2's C2 distance factor is wrong: u+v=2 eps r12, 2u+w=4 eps r2,
  2v+w=4 eps r1, hence C2=128 eps^5 r1 r2 r12, not 32 eps^5.
- The 33 and 177 counts require E=-eps^2; for fixed eps and general spectral E,
  they are 57 and 225. Keep this distinct from the fixed-eps physical-pencil theorem.
- “Through order 3” must distinguish unweighted exact output identities from
  admissible weighted approximation spaces; exact finite weighted jet inclusion is OPEN.
- Section 4 must require nested form-dense symmetry-correct spaces, physical
  dependency quotients, all nonhead weak rows, and a nonzero bounded head.
  Admissibility alone is insufficient. Residual zeros refer to nearby Ritz branches
  visible to the head; the inverse formula is understood through its singular limit.
- Abstract/Section 5: geometric decay is excluded only for the full plain array;
  raw tail experiments do not prove universal lack of dichotomy or Pincherle separation.
- “Algebraic convergence,” slopes -8 to -15, and order-3 improvement are numerical
  observations, not rate theorems. Use corrected r12 V data and compare basis sizes.
- Do not turn Morgan's convergence result into physical-eigenfunction representability,
  or FHOS regularity into analytic triple-point coefficient functions (ANSATZ.md).
