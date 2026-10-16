"""Data for the codex_psi30 figures: the Codex analytic series for psi30, by partial wave.

Run from the Codex workspace root (read only; nothing there is modified):
    cd ~/Documents/Codex/2026-09-05/how-much-compu
    PYTHONPATH=work/python-deps python3 ~/FockExpansion.jl/examples/sciviz/codex_series_data.py

Writes CSV files to examples/sciviz/output/codex_series/. The point-independent pieces
(Gauss coefficients g_e, multipole constants C_l, beta moments at t=1) are cached in
memory, so one evaluation at L=LMAX yields every truncation L <= LMAX at once.
"""
import functools, os, sys, time
from fractions import Fraction

ROOT = os.getcwd()
sys.path.insert(0, os.path.join(ROOT, "outputs/psi30_analytic"))
import mpmath as mp
import known_components as kc
import representation as rep
from full import psi30 as psi30_full

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "codex_series")
os.makedirs(OUT, exist_ok=True)
mp.mp.dps = 25
Z, E, A21 = 2, mp.mpf("-2.903724377034119598311159245194404446696925310"), mp.mpf("0.47674787900")

# Memoize the point-independent pieces (runtime only; the Codex files are untouched).
rep.gcoeff = functools.lru_cache(None)(rep.gcoeff)
rep.multipole = functools.lru_cache(None)(rep.multipole)
rep.low_constant = functools.lru_cache(None)(rep.low_constant)
kc._moment = functools.lru_cache(None)(kc._moment)


def terms(alpha, theta, L, N):
    """Elementary part and per-l terms of the three Legendre series, scaled as in full.py."""
    alpha, theta = mp.mpf(alpha), mp.mpf(theta)
    mu = mp.cos(theta)
    a = A21 - Z * (mp.pi + 4) / (9 * mp.pi)
    x = min(alpha, mp.pi - alpha)
    t = mp.tan(x / 2)
    base = kc.elementary(alpha, mu, Z, E, a) + Z**2 * kc.wave(1, t, "2b") * mu
    known = [Z**2 * (kc.wave(l, t, "2c") + kc.wave(l, t, "2e")) * mp.legendre(l, mu) for l in range(L + 1)]
    d2 = [Z**2 * rep.wave(l, x, "2d", N=N) * mp.legendre(l, mu) for l in range(L + 1)]
    # (1d) lives in the interelectron frame: polar angle y, even waves in mu_ee.
    cy = max(mp.mpf(-1), min(mp.mpf(1), -mp.sin(alpha) * mu))
    y, sy = mp.acos(cy), mp.sqrt(max(mp.mpf(0), 1 - cy * cy))
    muee = max(mp.mpf(-1), min(mp.mpf(1), mp.cos(alpha) / sy)) if sy else mp.mpf(0)
    panel, xe = ("A", y) if y <= mp.pi / 2 else ("B", mp.pi - y)
    d1 = [Z * rep.wave(l, xe, "1d", panel, N) * mp.legendre(l, muee) / mp.sqrt(2) if l % 2 == 0 else mp.mpf(0)
          for l in range(L + 1)]
    return base, known, d2, d1


def write(name, header, rows):
    with open(os.path.join(OUT, name), "w") as f:
        f.write(",".join(header) + "\n")
        for r in rows:
            f.write(",".join(str(v) for v in r) + "\n")
    print("wrote", name, len(rows), "rows", flush=True)


f = lambda v: mp.nstr(v, 20)

# Probe points (alpha, theta): generic, near e-e coalescence (alpha=pi/2, theta->0),
# near the nucleus (alpha->0), collinear opposite side (theta->pi).
PROBES = [("generic", 0.7, 0.9), ("equator", 1.5, 1.5), ("near_r12_0.30", mp.pi / 2, 0.30),
          ("near_r12_0.10", mp.pi / 2, 0.10), ("near_r12_0.03", mp.pi / 2, 0.03),
          ("near_r1_0.05", 0.05, 1.2), ("collinear", 1.0, mp.pi - 0.02)]
LMAX, NREF = 80, 240

# Sanity check against the unmodified Codex evaluator.
b, k, d2, d1 = terms(0.7, 0.9, 40, 200)
mine = b + sum(k) + sum(d2) + sum(d1)
ref = psi30_full(0.7, 0.9, Z=Z, E=E, a21=A21, L=40, N=200)
print("check vs full.psi30:", mp.nstr(mine - ref, 5), flush=True)

rows, timing = [], []
for name, al, th in PROBES:
    t0 = time.time()
    b, k, d2, d1 = terms(al, th, LMAX, NREF)
    timing.append((name, LMAX, NREF, time.time() - t0))
    rows.append((name, f(al), f(th), "base", 0, f(b)))
    for comp, arr in (("known", k), ("d2", d2), ("d1", d1)):
        rows += [(name, f(al), f(th), comp, l, f(v)) for l, v in enumerate(arr)]
    print(name, "psi30 =", mp.nstr(b + sum(k) + sum(d2) + sum(d1), 18), flush=True)
write("probe_terms.csv", ["point", "alpha", "theta", "component", "l", "term"], rows)

# Inner (Gauss) truncation N at fixed L=40; coefficient prefixes do not depend on N.
rows = []
for name, al, th in PROBES:
    for N in (10, 20, 30, 40, 60, 80, 100, 140, 200, 300):
        b, k, d2, d1 = terms(al, th, 40, N)
        rows.append((name, N, f(b + sum(k) + sum(d2) + sum(d1))))
write("probe_N.csv", ["point", "N", "psi30"], rows)

# Grid: partial sums S_L for every L <= LG, interior grid as in common.jl's interior(n).
NG, LG, NGRID = 41, 40, 200
eps = mp.mpf("1e-3")
nodes = [eps + (mp.pi - 2 * eps) * i / (NG - 1) for i in range(NG)]
rows, t0 = [], time.time()
for i, al in enumerate(nodes):
    for j, th in enumerate(nodes):
        b, k, d2, d1 = terms(al, th, LG, NGRID)
        s = b
        for L in range(LG + 1):
            s += k[L] + d2[L] + d1[L]
            rows.append((i, j, L, f(s)))
    print("grid row", i, "%.1fs" % (time.time() - t0), flush=True)
timing.append(("grid_per_point", LG, NGRID, (time.time() - t0) / NG**2))
write("grid_partial_sums.csv", ["ialpha", "itheta", "L", "psi30_L"], rows)

# Slow region: the line alpha = pi/2 (t = 1 for the nuclear-frame waves), to high L.
# Large l needs more digits (2^(l/2) cancellation) and a longer Gauss series.
mp.mp.dps, LL, NL = 40, 200, 360
rows = []
for th in ("0.03", "0.1", "0.3", "0.6", "1.0", "1.5", "2.2", "3.0"):
    b, k, d2, d1 = terms(mp.pi / 2, mp.mpf(th), LL, NL)
    rows.append((th, "base", 0, f(b)))
    for comp, arr in (("known", k), ("d2", d2), ("d1", d1)):
        rows += [(th, comp, l, f(v)) for l, v in enumerate(arr)]
    print("line theta", th, flush=True)
write("line_terms.csv", ["theta", "component", "l", "term"], rows)
mp.mp.dps = 25

# Warm cost (caches filled) per truncation L at a generic point.
for L in (5, 10, 20, 40, 80):
    terms(0.71, 0.9, L, 200)
    t0 = time.time()
    for k in range(5):
        terms(0.7 + 0.01 * k, 0.9, L, 200)
    timing.append(("warm", L, 200, (time.time() - t0) / 5))

# Cold cost of the unmodified evaluator (no caches) for the cost comparison.
import importlib, full, representation, known_components
for L, N in ((10, 60), (20, 140), (40, 200)):
    importlib.reload(known_components); importlib.reload(representation); importlib.reload(full)
    t0 = time.time()
    full.psi30(0.7, 0.9, Z=Z, E=E, a21=A21, L=L, N=N)
    timing.append(("full_py_uncached", L, N, time.time() - t0))
write("timing.csv", ["what", "L", "N", "seconds"], timing)

# Certified R1(l): exact rational start values and the step-two first difference.
def d_exact(l):
    c = Fraction(2 * l + 3, 2)
    g = Fraction(1)  # Gamma(c+2)/(sqrt(pi) Gamma(c+7/2)) = Gamma(l+7/2)/(sqrt(pi) Gamma(l+5))
    for j in range(l + 3):
        g *= Fraction(2 * j + 1, 2)
    for j in range(1, l + 5):
        g /= j
    return -4 * (8 * c * c + 18 * c + 1) * g / (3 * (2 * c - 1) * (2 * c + 1))


R1 = {2: Fraction(29, 192), 3: Fraction(1, 24)}
for l in range(2, 199):
    R1[l + 2] = R1[l] + d_exact(l)


def r1_quad(l):  # independent check: C1(c) = pi gamma(c)[log 2 + R1(l)] by quadrature
    c = mp.mpf(l) + mp.mpf(3) / 2
    P = lambda z: 1 - 6 * z / c + 12 * z * z / (c * (c + 1))
    F = lambda z: mp.hyp2f1(mp.mpf(5) / 2, mp.mpf(-3) / 2, c, z)
    C1 = mp.quad(lambda s: 2 * P(1 - s * s) * F(s * s), [0, 1])  # z = s^2
    gam = 3 * mp.gamma(c) / (4 * mp.sqrt(mp.pi) * mp.gamma(c - mp.mpf(5) / 2) * c * (c + 1))
    return C1 / (mp.pi * gam) - mp.log(2)


rows = []
for l in sorted(R1):
    q = f(r1_quad(l)) if l <= 60 else ""
    rows.append((l, R1[l].numerator if l <= 12 else "", R1[l].denominator if l <= 12 else "",
                 f(mp.mpf(R1[l].numerator) / R1[l].denominator), f(mp.mpf(d_exact(l).numerator) / d_exact(l).denominator), q))
write("r1.csv", ["l", "num", "den", "R1", "d_l", "R1_quadrature"], rows)
