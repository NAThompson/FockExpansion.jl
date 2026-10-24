# Claude working notes

## Literature (checked 2026-10-21 by WebFetch of arXiv abstracts/ar5iv)
- FHOS 2009 CMP 289:291 (arXiv:0806.1004): Thm 1.1 is LOCAL and BINARY: psi = psi1 + |x1| psi2, psi_i real analytic,
  near one-electron–nucleus coalescence (others separated); similar near two-electron coalescence. Triple point NOT treated.
- FHOS 2005 CMP 255:183 (math-ph/0312060), "Sharp regularity results for Coulombic many-electron wave functions":
  psi = e^{F2+F3} phi, phi in C^{1,1}; F2 = Kato cusp terms, F3 = C0 sum Z (x_i.x_j) ln(|x_i|^2+|x_j|^2),
  C0 = (2-pi)/(12 pi) in units H = -Delta + V (x = 2 x_au). Converted: equals Fock psi21 R^2 ln R. Only the first log;
  higher logs NOT absorbed; regularity only C^{1,1}, not analyticity.
- So the hint "Psi = F1 + R^2 ln R F2 with F1,F2 analytic" is NOT a theorem as stated (to my reading). Morgan 1986:
  Fock expansion converges for small R, with ln^k R, k <= floor(n/2).
