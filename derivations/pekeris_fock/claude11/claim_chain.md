# claim_chain (Claude, 2026-10-24) — numerical, n=80/N=60/256-bit vs n=64/N=48
R1 side -1: no logs at X0 = -3pi/2 through psi40 (all A_m, m>=1, below 1e-29 at 128 bits; exactly dropped at 256 bits).
R2 side +1: log degree 1 at Xb = 3pi/2 for psi20 (Z^1), psi30 (Z^1, Z^2), psi40 (Z^1, Z^2, Z^3), psi41 (Z^2). Continuation agrees
   with the physical spectral interpolant just beyond the edge to 1e-49; far-edge log coefficients converge n64 vs n80 to 4e-37.
R3 chain: step-1 (N_a N_b) nonzero only for psi40 Z^2: B2(alpha) = 1/9 - (7/45) sin(2 alpha) (fit residual 5e-52), N_b N_a N_b = 0.
Files: extend.jl, chain.jl, chain_80.log, chain_64.log, analyze_chain.jl, chain_sectors.jls, b2fit.jl, guess_ode.jl.
