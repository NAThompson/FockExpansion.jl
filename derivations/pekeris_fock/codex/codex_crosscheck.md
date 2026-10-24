# Independent derivation and cross-check — 2026-10-21

## Provenance and reproducibility

`codex_derive.py`, `codex_operator.json`, `codex_recurrence.json`, and the passing derivation log were produced **before** opening any file under `sym/` or `julia/`. Only then were `sym/recurrence.json`, `sym/recurrence.py`, and `sym/verify_recurrence.py` read. No file under `julia/` was read. The initial independent algorithm uses the exact derivative-tail formula, rather than Claude's elimination of bare derivatives with the Laguerre differential equation.

Run from the workstream root:

    venv/bin/python -B codex/codex_derive.py > codex/codex_derive.log 2>&1
    venv/bin/python -B codex/codex_compare.py
    venv/bin/python -B sym/verify_recurrence.py quick

The last command is read-only with respect to Claude's files; `-B` prevents bytecode writes there. Its captured output is `codex_claude_verify.log`. No repairs to Claude's code were made.

## Independent derivation

Write a=r1,b=r2,c=r12 and d=(d_a,d_b,d_c). The S-state distance operator is

    H = -(1/2) [ d_a^2+d_b^2+2d_c^2
          + (2/a)d_a+(2/b)d_b+(4/c)d_c
          + (a^2+c^2-b^2)/(a c) d_a d_c
          + (b^2+c^2-a^2)/(b c) d_b d_c ]
        - Z/a - Z/b + lambda/c.

This follows by applying the Cartesian Laplacians to functions of the three distances; the mixed terms are the dot products of their distance gradients. Code constructs the symmetric second-order matrix G, where off-diagonal mixed coefficients are halved, and drift h=(2/a,2/b,4/c).

For x=(u,v,w), the constant Jacobian J=dx/d(a,b,c) is

    eps * [ -1  1  1
             1 -1  1
             2  2 -2 ].

The inverse is a=(2v+w)/(4eps), b=(2u+w)/(4eps), c=(u+v)/(2eps). Let Q=(u+v)(2u+w)(2v+w)=32eps^3 abc. Our normalization is **exactly** COLLAB.md's P=-Q exp(sigma/2)(H-E)exp(-sigma/2), including its outer eps factor. There is no residual 2eps normalization difference.

Before conjugation, second-order matrix and drift are Q JGJ^T/2 and Q Jh/2, and zeroth term is Q(Z/a+Z/b-lambda/c+E). Substitution makes every coefficient polynomial. Conjugation replaces each d_i by d_i-1/2. The saved differential operator uses one coefficient for each derivative multi-index, so mixed entries of M are added twice.

For any symmetric M, drift b, and sigma=u+v+w, the product rule gives

    C1 = 2sigma sum_(i,j) M_ij d_j + sigma sum_i b_i - sum_(i,j) M_ij,
    C2 = sum_(i,j) M_ij = 4eps^2 Q.

Thus the jth log-sector equation is sigma^2 P F_j +(j+1) C1 F_(j+1) +(j+1)(j+2) C2 F_(j+2)=0.

## Exact Laguerre conversion and tests

The independent conversion starts from

    x L_k = -(k+1)L_(k+1)+(2k+1)L_k-k L_(k-1),
    D^b L_q = (-1)^b sum_(k=0)^(q-b) binomial(q-k-1,b-1) L_k, b>=1.

For x^a D^b, repeated multiplication gives a finite set of intermediate shifts. The derivative formula gives their exact equation-index contribution for each offset. Terms are combined before deciding which offsets vanish. All derivative orders are <=2; beyond the maximum coordinate degree the potential infinite tails are polynomials of degree <=1 in an offset. The code checks two or more consecutive zero slices in that region, with all other offset categories included; this certifies cancellation, rather than simply assuming a finite stencil. It checks offsets through +8; the largest coefficient coordinate degree is 5. At nonnegative equation indices the multiplication factors kill paths from negative intermediate Laguerre indices; negative source coefficients are defined as zero.

For each of P, sigma^2 P, C1, C2, direct symbolic differentiation on the three basis sources (0,0,0), (1,2,0), (3,1,2) is compared with the entire predicted output polynomial. Projection uses the independent integral identity integral exp(-x) x^p L_i(x) dx=(-1)^i p! binomial(p,i). All **12 exact polynomial-action tests pass**, using eps=3/2,E=-7/3,Z=2,lambda=4/5 so the test does not assume E=-eps^2 and does not silently set lambda=1. Zero-index boundaries are included. Three exact symbolic product-rule tests apply P to log(sigma)^k(1+uv+w^2), k=1,2,3, and compare to the sector identity; all pass.

The JSON coefficients are in equation-index convention:

    at (l,m,n): sum_d coefficient_d(l,m,n) A(l+dl,m+dm,n+dn)=0.

`codex_recurrence.json` stores **general E** for all four keys `P`, `S2P`, `C1`, `C2`. Substitute E=-eps**2 to obtain the Pekeris specialization. Symbols are `eps`, `E`, `Z`, `lambda`, `l`, `m`, `n`. Because `lambda` is a Python keyword, a consumer using SymPy's Python-expression parser must alias that token before parsing; `codex_compare.py` demonstrates a token-bounded replacement by `lam`. This is a serialization issue, not a change in the coupling symbol.

| operator | general offsets | shell range | E=-eps^2 offsets | shell range |
|---|---:|---|---:|---|
| P | 57 | [-3,3] | 33 | [-2,2] |
| sigma^2 P | 225 | [-5,5] | 177 | [-4,4] |
| C1 | 117 | [-4,4] | 117 | [-4,4] |
| C2 | 57 | [-3,3] | 57 | [-3,3] |

## Term-by-term comparison with Claude

Reference `sym/recurrence.json` SHA-256 at comparison:

    953df573062bce6c24ea2cc4948967bba131a3132642cdbd9335d3190ad7324d

Claude's `sigma2P` is our `S2P`, and `epsilon` is our `eps`. His JSON stores fixed-energy sectors, plus `P_general`; therefore the comparison specializes our sector coefficients before comparing. Symbols are parsed into the same SymPy objects to remove assumption/name differences. The exact expanded difference is checked for **every offset in the union**, with a missing coefficient interpreted as zero.

| comparison | coefficients checked | nonzero symbolic differences |
|---|---:|---:|
| P general | 57 | 0 |
| P fixed | 33 | 0 |
| sigma^2 P fixed | 177 | 0 |
| C1 | 117 | 0 |
| C2 | 57 | 0 |
| total | **441** | **0** |

Machine-readable results and empty discrepancy maps are in `codex_comparison.json`. In particular **there is no independent C1 or C2 disagreement**.

## Reported sector-machine-test issue

The current `sym/verify_recurrence.py quick` gives:

    seed 0, eps^2=-E: False, sector 0: 209 Laguerre coefficients; support mismatch 0 value mismatch 0
    seed 0, eps^2=-E: False, sector 1: 209 Laguerre coefficients; support mismatch 0 value mismatch 0
    seed 0, eps^2=-E: False, sector 2: 212 Laguerre coefficients; support mismatch 0 value mismatch 0
    ALL OK

Thus the mismatch described in TASK_codex_1.md is **not reproducible with the currently present files**. No cause for the earlier failure is asserted and no fix is claimed. The independent tests and coefficient comparison support the current stencils, including both log couplings. This run is the quick seed-0 case; it is not a claim that every optional seed in Claude's script was run.

An additional exact algebra check found zero discrepancies in the formal symmetry relation P_d(i)=P_(-d)(i+d) for all 57 general-energy offsets. This supports the generalized symmetric operator-pencil interpretation discussed in NOTES_codex.md; it is not a self-adjointness/domain proof.
