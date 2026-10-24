# Small serial verification — 2026-10-22

Command: `PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase5/goal1_verify.py`

```text
input_sha256 928d8ca17b031c83d86902c4cece68d639ae281b5361d7a9787ec59ba5458564
precision 45 cutoffs 2,4,6; largest dimension 50; serial
N 2 z -5.0 Schur_residual 1.0907541904834e-48 direct_difference 6.2542694216157e-49 normalized_source 32.053687898669
N 2 z (-2.9037243770341 + 0.05j) Schur_residual 1.2524544595011e-48 direct_difference 4.2099887019956e-47 normalized_source (0.13873644063958 - 1.8958995019376j)
N 2 size 7 E0_tail_positive definite; defect 0.017337815181131 raw_P_dropped_row -0.017337815181131
N 4 z -5.0 Schur_residual 4.6928130511269e-46 direct_difference 1.0477237663561e-47 normalized_source 28.271085035311
N 4 z (-2.9037243770341 + 0.05j) Schur_residual 9.3013813247754e-46 direct_difference 7.1214657108419e-48 normalized_source (0.13563334481771 - 1.9000014248675j)
N 4 size 22 E0_tail_positive definite; defect 0.0013538011999736 raw_P_dropped_row -0.0013538011999736
N 6 z -5.0 Schur_residual 1.1249163440768e-45 direct_difference 3.2019721226905e-48 normalized_source 27.059594148238
N 6 z (-2.9037243770341 + 0.05j) Schur_residual 7.9813257315475e-46 direct_difference 7.7146582993162e-47 normalized_source (0.13982833921063 - 1.9004513494993j)
N 6 size 50 E0_tail_positive definite; defect 0.00018428754978695 raw_P_dropped_row -0.00018428754978695
PASS: Hermiticity, Gram/shifted/tail positivity, bandwidth, backward Schur,
normalized-tail equivalence, all retained tail equations, decreasing E0 defect.
Finite checks support algebra only; the infinite-limit proof is in goal1_theorem.md.
peak_RSS_MiB 65.97
```
