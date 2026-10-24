# Pekeris recurrence with exact Fock structure

This directory holds a research record from a Claude–Codex collaboration (October 2026). It shows that the helium ground state is the selected limit of an explicit finite recurrence. The recurrence has the exact Fock structure built in through order 3.

Start with:

- `THEORY.md`: the current theorems, the verified claims and the open problems.
- `COLLAB.md`: the decision log, with who did what and each co-sign or dissent.
- `ANSATZ.md`: the global ansatz e^{-σ/2} Σ_k ln^k σ · Σ_i g_i F_{k,i}, with σ = u+v+w.

Each phase has two folders, `claudeN/` (Claude's derivations and numerics) and `codex/phaseN/` (Codex's independent work):

| Phase | Content |
|---|---|
| 1 (`sym/`, `julia/`) | Pekeris recurrence plus the ln σ sectors, and its verification |
| 2–3 | Order-2 Fock atoms (minimal rank-5 set) and ψ30 numerics |
| 4 | Fock-constrained basis and the a21 readout |
| 5–6 | Zero-tail backward theorem; the order-3 degree blow-up in the cyclic basis |
| 7–8 | Structured order-3 basis: the classical block and four elliptic panels, giving 2,577,274 stencil terms |

The compiled stencil generator is `claude8/stencil3.jl`. It regenerates the order-3 stencils from `claude8/*_terms.txt` and the row-factor files; the expanded stencils, about 6 GB per panel, are not stored here.

`prior/` holds the earlier boundary-matching studies:

- the exterior inverse-factorial proof;
- the Steklov–Poincaré a21 extraction.

Large caches (`*.jls`, `*.pkl`, expanded stencils) and third-party literature copies are left out.
