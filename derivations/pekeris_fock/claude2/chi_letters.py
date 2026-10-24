"""chi (Z^1 part of psi20, repo src/FockExpansion.jl) rewritten in 'letters' of the field K = Q(r1,r2,r12)(R,P),
R = sqrt(r1^2+r2^2) = hyperradius, P = sqrt(2R^2 - r12^2) = |x1+x2|; signed d = r1-r2, s_h = r1+r2.
Weight-1 atoms: alpha = 2 atan(r2/r1), beta = asin(v), v = (R^2-r12^2)/R^2, and logs of letters.
Weight-2: DL = Cl2(alpha-beta)+Cl2(pi+beta-alpha)-Cl2(pi-alpha-beta)-Cl2(alpha+beta)."""
import mpmath as mp

mp.mp.dps = 60
cl2 = lambda t: mp.clsin(2, t)


def chi_letters(r1, r2, r12):
    r1, r2, r12 = map(mp.mpf, (r1, r2, r12))
    R = mp.sqrt(r1**2 + r2**2); P = mp.sqrt(2 * R**2 - r12**2)
    sh, d = r1 + r2, r1 - r2
    s = 2 * r1 * r2 / R**2; c = d * sh / R**2; v = (R**2 - r12**2) / R**2
    g = r12 * P / R**2; sig = sh / R; xi = r12 / R
    al = 2 * mp.atan(r2 / r1); be = mp.asin(v)
    h = -(3 * mp.pi + 10 - 16 * cl2(mp.pi / 2)) / (24 * mp.pi)
    L5 = r12 * P + sh * d
    lnsx = mp.log(sh + r12) - mp.log(R)                      # log(sigma + xi)
    Lq = 2 * mp.log(d + r12) - mp.log(L5)                   # log((c+sigma xi)^2/(sigma^2 (g+c)))
    logs = (be * (2 * mp.log(2 * r1 * r2) - 2 * mp.log(L5))
            + (be + al) * (mp.log(sh - r12) + mp.log(sh + r12) - mp.log(2 * r1 * r2))
            + (be - al) * (mp.log(r12 - d) + mp.log(r12 + d) - mp.log(2 * r1 * r2)))
    DL = cl2(al - be) + cl2(mp.pi + be - al) - cl2(mp.pi - al - be) - cl2(al + be)
    return ((g + s - 2 * sig * xi) / 6 + h * v
            + (-2 * mp.pi * v * lnsx + mp.pi * c * Lq + 2 * g * be + c * (logs - DL)) / (6 * mp.pi))


if __name__ == '__main__':
    for line in open('chi_ref.txt'):
        a, b, c, ref = line.split()
        val = chi_letters(float(a), float(b), float(c))  # Julia used the Float64 inputs
        print(a, b, c, 'diff =', mp.nstr(val - mp.mpf(ref), 5))
