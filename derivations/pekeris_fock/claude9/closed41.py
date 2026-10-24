# Assemble the exact form of R^4 psi41^(2) from the identified coefficients (module {1, alpha, L, U, V'}).
import re, sympy as S
r1, r2, r12, pi, R = S.symbols('r1 r2 r12 pi R', positive=True)
a_h, b_h = S.symbols('a_h b_h')   # harmonic constants (numerical)
terms = {}
for line in open('coef5_ident.txt'):
    lab, val, rel = [s.strip() for s in line.split('|')]
    if rel == 'NONE':
        continue
    c = [int(t) for t in re.findall(r'-?\d+', rel)]
    x = -(c[1] + c[2] * pi + c[3] * pi**2 + c[4] / pi + c[5] / pi**2) / c[0]
    atom, mono = lab.split(' ', 1)
    i, j, k = [int(t) for t in re.findall(r'\^(\d+)', mono)]
    m = i + j + k
    terms[atom] = terms.get(atom, 0) + S.nsimplify(x) * r1**i * r2**j * r12**k * R**(4 - m - (1 if atom == 'V' else 0))
# even polynomial part: rest + a Y40 + b Y42 (R^4-homogenized)
rest = (5 - 3 / pi - 14 / pi**2) / 360
Y40 = 4 * (r1**2 - r2**2)**2 - R**4
Y42 = S.Rational(3, 2) * (R**2 - r12**2)**2 - 2 * r1**2 * r2**2
terms['1'] = terms.get('1', 0) + rest * R**4 + a_h * Y40 + b_h * Y42
for k, v in terms.items():
    v = S.factor(S.expand(v.subs(R, S.sqrt(r1**2 + r2**2))))
    print(f'{k}: {v}')
    terms[k] = v
import pickle; pickle.dump(terms, open('closed41.pkl', 'wb'))
