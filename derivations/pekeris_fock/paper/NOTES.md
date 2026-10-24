# NOTES (drafting decisions)

Codex brevity memo: BREVITY_codex.md. Followed except where noted.

Followed:
- Dropped the cyclic-vector reconstruction story to one sentence (the degree >= 351 artefact; the user asked that it be stated).
- Dropped the offsets table for the plain sector operators (counts go in a sentence).
- One construction table as the central object; no figure.
- Separated the two results (explicit recurrence vs. weak-pencil selection theorem); their identification is stated as OPEN.
- "Built-in decay" reworded: the exponential is a gauge; decay is imposed by form-domain admissibility.
- C2 = 128 eps^5 r1 r2 r12 (corrected factor).
- 33/177 counts tied to E = -eps^2; 57/225 at general E.
- Supplement material referred to the repository paths instead of a separate supplement file.

Disagreements (reasons):
- Codex suggested dropping the energy table and the a21 estimate. The user explicitly asked for the E-error table and a21 as a remark.
  Kept a compact table (corrected r12*V data, with basis sizes, labelled numerical) and one sentence on a21 with its empirical scale.
- Codex suggested reducing "no exponential dichotomy" to one sentence. The user asked for the Jaffe/Leaver comparison; kept one short
  paragraph (no plot), with the proved part (no geometric decay of the full plain array) separated from the numerical observation.
- Open problems: kept a short list (5 items) rather than only two, because the user asked for open problems; each is one line.

## After Codex review (REVIEW_codex.md)
Applied all MUST-FIX items 1-18, 20-21: "propose the ansatz" wording; count scoped to five unweighted blocks at E=-eps^2, Z=2;
abstract/numerics reworded to "consistent with algebraic"; finite-tower statement with all-sector domain hypothesis; sector equations
used formally; C2 "proportional to" the density; admissibility stated for exponentially damped products; cyclic degree evidence stated as a
finite-field reconstruction failure through total degree 350; involution with its Jacobian; singlet-S scope of the core; local scope of
the Ritz-zero statement and the Schur/min-max sketch; "spectral hypothesis" instead of "cited input"; |g_N|; corrected entire-function
implication; Pincherle definition; energy-table caption describes the actual bases; a21 convention stated; repository URL + commit.
Applied optional cuts: timings/monomial totals/storage, duplicate Miller-vs-Ritz numbers, k^-4 sentence, duplicate density sentence.
Not applied: item 22 (credit Liverts-Barnea etc. for the order-2/3 coefficient formulas) to save space; the convention is named in the
a21 remark. Item 23: Codex confirmed Pekeris 1958/1959 and Koutschan-Zeilberger via APS and the authors' page; TODOs removed for those.
