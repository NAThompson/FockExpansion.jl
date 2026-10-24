# OUTLINE (candidate content, before brevity cuts)

Working title: "Helium's ground state as the selected limit of a finite recurrence with exact Fock structure"
Target: 4-6 pages, letter style (revtex4-2 or article two-column-ish). Author: Nick Thompson (paper/ in repo has an empty \author{},
so fall back to the default), acknowledgement of Claude (Anthropic) and Codex (OpenAI).

## 0. Abstract (~120 words)
Global ansatz e^{-sigma/2} sum_k ln^k sigma sum_i g_i F_{k,i}; finite recurrences through order 3 (2,577,274 offsets, +-16 shells);
zero-tail backward selection theorem; no exponential dichotomy; algebraic convergence. Honesty labels.

## 1. Introduction (~0.6 p)
- Special-function wish list: (i) exact Fock structure at the triple point, (ii) built-in decay, (iii) explicit finite recurrence,
  (iv) physical solution selected by a recurrence criterion.
- Pekeris 1958/59: 33-term recurrence in perimetric coords, but analytic F cannot carry Fock's R^2 ln R (Fock 1954; FHOS 2005;
  Frankowski-Pekeris 1966 added logs variationally; Morgan 1986 convergence of Fock series).
- Contributions list (4 bullets) with labels.

## 2. Ansatz and log sectors (~0.7 p)
- Perimetric u,v,w; sigma = 2 eps (r1+r2); P operator; sector equations sigma^2 P F_j + (j+1) C1 F_{j+1} + (j+1)(j+2) C2 F_{j+2} = 0.
- C2 = 4 eps^2 (u+v)(2u+w)(2v+w) = 32 eps^5 r1 r2 r12 (check factor).
- Table 1: offsets P 33, sigma^2 P 177, C1 117, C2 57; shell ranges. VERIFIED (two independent derivations, 441 coeffs).
- Log-free class restriction, gauge, finite-K tower has no log content (PROVED under domain hypothesis).

## 3. Atoms (~1 p)
- Order 2: chi; minimal rank-5 module {1, L, alpha, V, U} (PROVED, monodromy). r12 V for form domain. Finite recurrence,
  half-bandwidth 10; 25,955 offsets in the 8-atom presentation.
- Order 3: cyclic vector of {1,chi,R^3 psi30}: rank 45 candidate, line degree >= 351 (artefact). Structured basis: classical block
  (18-state, 21 sectors) + four elliptic panels (rational 55-state -> 23 reachable states over Q(pi)(r1,r2,r12)), fixed-endpoint
  globality (involution s -> h/s). Table 2: per-block sectors/offsets/bandwidth: classical 21/215,406/12; panels 1,2 23/562,570;
  3,4 23/618,364; total 2,577,274; +-16; row factors deg <= 12. Direct-operator checks 1e-30..1e-54; Codex certified selected blocks.
- Admissibility: classical unweighted (PROVED), panels with delta^4 (PROVED sufficient), powers <1 fail, minimal power in {1..4} OPEN.

## 4. Selection theorem (~1 p)
- Theorem (zero-tail backward): A000 bounded on form domain, coercivity on ker A000, Cea; dropped row g_N = 1/[K_N^{-1}]_00 -> 0,
  zeros = Ritz energies. Proof sketch: density (holomorphic tube / Fourier), Hardy, gap.
- Extension to nested weighted structured spaces (PROVED under same cited input).
- Not a formal cleared-recurrence theorem (OPEN), no N-rate.

## 5. No exponential dichotomy (~0.4 p)
- Numerics: tail data propagate undamped (Fig. goal1_separation). Leaver/Jaffe: double root 1 even classically; Pincherle needs no
  geometric dichotomy. PROVED: geometric decay impossible due to R^2 ln R. No Pincherle classification.

## 6. Numerics (~0.5 p)
- Table 3: E - E0 vs N for polynomial, +log, +log^2, min5 (r12 V), min5+psi30. Algebraic (log-log slope ~ -8 to -15, semilog convex).
- Figure: convergence3.png (or goal1_separation.png).
- Remark: a21 ~ 0.47674788 (empirical scale 1e-8).

## 7. Open problems (~0.3 p)
- N-rate; formal cleared-recurrence version; classification of formal solutions; minimal panel weight; unweighted jet inclusion;
  minimal ranks; order 4 (psi40, psi41).

## Acknowledgements; references (~0.4 p)
Fock 1954; Pekeris 1958, 1959; Frankowski-Pekeris 1966; FHOS 2005 (and 2009?); Morgan 1986; Leaver 1985, 1986; Koutschan-Zeilberger 2011;
Klahn-Bingel 1977.

## Supplement candidates
prior/ boundary matching; all Codex/Claude cross-check logs; per-block count tables; atom formulas for chi, U.
