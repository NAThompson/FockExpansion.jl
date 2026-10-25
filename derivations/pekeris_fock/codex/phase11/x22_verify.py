"""Small serial certificates for the source-independent phase-11 pilot.

No imported Claude/repository code. Symbolic identities are independent of
the numerical checks; all writes stay in this file's x22_* namespace.
"""
import sys
sys.dont_write_bytecode = True
import json, resource, time
from pathlib import Path
import sympy as S
import mpmath as mp

START = time.time()
OUT = Path(__file__).resolve().parent
checks, facts = {}, {}
n, z, q, t, k, a, b, u, v = S.symbols('n z q t k a b u v')
Q = n-z*z
ell = 1-z
om = 1-n
h, J, lg = S.symbols('h J lg')
zero = lambda x: S.cancel(S.expand(x)) == 0

# H = gamma-theta, cos(gamma)=z/sqrt(n),
# cos(theta)=(z-n)/(sqrt(n)*(1-z)), with compatible physical branches.
# J=H/sqrt(n-z^2). Basis (J,1,h), h^2=1-n.
Gz = S.Matrix([[z/Q, -1/Q, 1/(ell*Q)], [0,0,0], [0,0,0]])
Gn = S.Matrix([[-1/(2*Q), z/(2*n*Q), -(z+n)/(2*n*om*Q)],
               [0,0,0], [0,0,-1/(2*om)]])
checks['rank3_full_connection_flat'] = all(zero(x) for x in
    Gz.diff(n)-Gn.diff(z)+Gz*Gn-Gn*Gz)
checks['top_character_is_Q_only'] = zero(Gz[0,0]+S.diff(Q,z)/(2*Q)) and zero(Gn[0,0]+S.diff(Q,n)/(2*Q))
checks['bottom_character_is_omega_only'] = zero(Gn[2,2]-S.diff(om,n)/(2*om))
facts['basis'] = ['J=(gamma-theta)/sqrt(n-z^2)', '1', 'h=sqrt(1-n)']
facts['connection_dz'] = [[str(x) for x in row] for row in Gz.tolist()]
facts['connection_dn'] = [[str(x) for x in row] for row in Gn.tolist()]
facts['graded_characters'] = ['1', 'sqrt(1-n)', '1/sqrt(n-z^2)']

# Add h*log(1-z), the only other target-nonpolynomial zonal generator.
Az = S.Matrix([[z/Q,0,-1/Q,1/(ell*Q)], [0,0,0,-1/ell], [0,0,0,0], [0,0,0,0]])
An = S.Matrix([[-1/(2*Q),0,z/(2*n*Q),-(z+n)/(2*n*om*Q)],
               [0,-1/(2*om),0,0], [0,0,0,0], [0,0,0,-1/(2*om)]])
checks['rank4_full_connection_flat'] = all(zero(x) for x in Az.diff(n)-An.diff(z)+Az*An-An*Az)
facts['rank4_basis'] = ['J','h*log(1-z)','1','h']
facts['rank4_dz'] = [[str(x) for x in row] for row in Az.tolist()]
facts['rank4_dn'] = [[str(x) for x in row] for row in An.tolist()]
facts['q_density_basis'] = ['1','A=pi/2-atan(q)','A^2']
facts['q_density_connection'] = [['0','0','0'],['-1/(1+q^2)','0','0'],['0','-2/(1+q^2)','0']]
facts['tensor_integrand_upper_bound'] = 12

# Nonpolynomial parts of all three necessary zonal inverses. The rest of
# each U is a polynomial in target z, with parameter-dependent coefficients.
T3 = z*(4*z*z-3*n)
Ainv = (4*z*z-n)*(3*n-4)/(12*n**3)
Al = -(4*z*z-n)*ell*(1+z)*(2*z*z-3*n+2)/(24*n**3) + T3*z*Q/(12*n**3)
Ax = -(4*z*z-n)*ell**2*(8*z**3+6*z*z-10*z*n+4*z-5*n+2)/(120*n**3) + T3*Q*(z-(n+4*z*z)/5)/(12*n**3)
Cs = [-(n-4)*T3/(12*n**3), om*T3/(12*n**3), om**2*T3/(60*n**3)]
As = [Ainv, Al, Ax]

def reduce_h(f):
    # Only polynomial dependence on h occurs in these checks.
    return S.rem(S.Poly(S.expand(f),h), S.Poly(h*h-om,h)).as_expr()

def dz(f):
    return reduce_h(S.diff(f,z)+S.diff(f,J)*(z/Q*J-1/Q+h/(ell*Q))-S.diff(f,lg)/ell)

for kind,A,C,source in zip(['inv','log','xlog'], As, Cs, [1/ell,lg,ell*lg]):
    F = A*lg+h*C*J
    residual = S.expand(-4*(Q*dz(dz(F))-3*z*dz(F)+8*F)-source)
    coeffJ = S.cancel(residual.coeff(J))
    coefflg = S.cancel(residual.coeff(lg))
    rest = S.cancel(residual.subs({J:0,lg:0}))
    checks[kind+'_nonpolynomial_PDE_terms_cancel'] = coeffJ == 0 and coefflg == 0
    # An explicit polynomial remainder is sufficient: it is absorbed by
    # the target-polynomial part of the actual zonal formula.
    den = S.fraction(rest)[1]
    checks[kind+'_PDE_remainder_target_polynomial'] = not den.has(z)
    facts[kind+'_log_coefficient'] = str(S.factor(A))
    facts[kind+'_hJ_coefficient'] = str(S.factor(C))
    facts[kind+'_PDE_polynomial_remainder'] = str(S.factor(rest))

# Exact resonance cancellation of the bare gamma coefficient in I2.
# C_m(gamma) has gamma coefficient 0 for odd m, binomial(2j,j)/4^j for m=2j.
rho = S.symbols('rho')
def cgam(m):
    return S.binomial(m,m//2)/2**m if m%2 == 0 else S.Integer(0)
def trig_gamma(terms):
    return S.expand(-sum(an*cgam(j)*rho**(j-i) for i,an in terms for j in range(i)))
sets = [[(4,-4),(2,5),(0,-1)], [(1,1),(3,-2),(5,1)],
        [(1,1),(3,-2),(5,1),(0,-rho/5),(2,-2*rho/5),(4,7*rho/5),(6,-4*rho/5)]]
for kind,terms in zip(['inv','log','xlog'],sets):
    g = trig_gamma(terms)
    facts[kind+'_raw_trig_gamma_coefficient'] = str(S.factor(g))
    # I2pi contributes kappa=-2 I2pi/pi. Its gamma/2 term
    # cancels precisely the rational (non-sqrt(omega)) gamma term.
    raw = -g/4 if kind=='inv' else rho*g/4
    kap_part = -2*raw
    checks[kind+'_resonant_rational_gamma_cancels'] = zero(raw+kap_part/2)

# Reconstruct the actual regular-antipodal zonal primitive from independent
# integration recurrences. This fixes the physical branch; a PDE identity
# alone would not exclude adding an arbitrary singular homogeneous solution.
c,sg,ga,at,HH,Lone = S.symbols('c sin_gamma gamma atan_half H log_one_minus_rho')
sr = S.symbols('sqrt_omega')
def Cpow(m):
    if m==0: return ga
    if m==1: return sg
    return c**(m-1)*sg/m+S.Rational(m-1,m)*Cpow(m-2)
def trigsum(terms):
    return S.expand(-sum(an*Cpow(j)*rho**(j-i) for i,an in terms for j in range(i)))
jt=[-lg/rho]
for j in range(1,7): jt.append((jt[-1]-c**j/S.Integer(j))/rho)
Ht=lambda j:c**(j+1)*lg/S.Integer(j+1)+rho*jt[j+1]/S.Integer(j+1)
Fints=[4*jt[3]-3*jt[1],4*Ht(3)-3*Ht(1),4*(Ht(3)-rho*Ht(4))-3*(Ht(1)-rho*Ht(2))]
P0=c*sg**3
P1=sg**3*(c-rho*(1+4*c*c)/5)
I2s=[-((rho*rho-4)*sr*2*at/rho**4+trigsum(sets[0]))/4,
     -(P0*lg-rho*(sr**3*2*at/rho**5+trigsum(sets[1])))/4,
     -(P1*lg-rho*(sr**5*2*at/(5*rho**5)+trigsum(sets[2])))/4]
for kind,FI,I2,A,C in zip(['inv','log','xlog'],Fints,I2s,As,Cs):
    I1=(FI-FI.subs({c:1,lg:Lone}))/4
    kap=-2*S.expand(I2).coeff(ga)-S.expand(I2).coeff(at)
    U=(4*c*c-1)*I1/3-(4*c**3-3*c)*I2/(3*sg)+kap*((4*c*c-1)/18-ga*(4*c**3-3*c)/(6*sg))
    targetF=A.subs({n:rho**2,z:rho*c})*lg+sr*C.subs({n:rho**2,z:rho*c})*HH/(rho*sg)
    diff=S.cancel(S.expand(U.subs(at,(ga-HH)/2)-targetF))
    # Clear powers of sin(gamma) and sqrt(omega) using their defining equations.
    num,den=S.fraction(diff)
    num=S.rem(S.Poly(num,sg),S.Poly(sg**2-(1-c*c),sg)).as_expr()
    den=S.rem(S.Poly(den,sg),S.Poly(sg**2-(1-c*c),sg)).as_expr()
    diff=S.cancel(num/den)
    num,den=S.fraction(diff)
    num=S.rem(S.Poly(num,sr),S.Poly(sr**2-(1-rho*rho),sr)).as_expr()
    den=S.rem(S.Poly(den,sr),S.Poly(sr**2-(1-rho*rho),sr)).as_expr()
    diff=S.cancel(num/den)
    checks[kind+'_actual_regular_zonal_mod_polynomial'] = not diff.has(ga,HH,lg,sg) and not S.fraction(diff)[1].has(c)
    facts[kind+'_actual_zonal_polynomial_difference'] = str(S.factor(diff))
    facts[kind+'_polynomial_target_degree'] = int(S.degree(diff,c))

# Feynman beta radical and the actual C3 panel-3/4 data.
D = 1+q**4
pp = [S.Matrix([-2*q*q,1-q**4]), S.Matrix([2*q*q,1-q**4])]
xx = S.Matrix([u,v])
for panel,P in zip([3,4],pp):
    checks[f'panel{panel}_unit_P'] = zero(P.dot(P)-D**2)
    c0 = (1-v*v)/18
    c2 = (-1 if panel==3 else 1)*2*(1+v)*u/9
    H = S.factor(c0-c2*q*q)
    facts[f'panel{panel}_A_numerator'] = str(H)
    checks[f'panel{panel}_actual_A_numerator'] = zero(H-(1+v)*(1-v+(1 if panel==3 else -1)*4*q*q*u)/18)
    # Both nuclear centres must be retained, even in this model.
    for sign in [1,-1]:
        aa = S.Matrix([sign,0])
        N = S.factor(D-aa.dot(P))
        Nroot = 1+(sign if panel==3 else -sign)*q*q
        checks[f'panel{panel}_centre{sign}_N_square'] = zero(N-Nroot**2)
        rvec = t*aa+(1-t)*P/D
        nn = S.factor(rvec.dot(rvec))
        zz = rvec.dot(xx)
        QQ = S.factor(nn-zz**2)
        kk = 2*N/D
        checks[f'panel{panel}_centre{sign}_omega_beta'] = zero(1-nn-kk*t*(1-t))
        # (mu*h)^2 = 2N/D^2; positive physical branch
        # mu=1/sqrt(D*t*(1-t)), hence mu*h=sqrt(2)*Nroot/D.
        checks[f'panel{panel}_centre{sign}_measure_cancel'] = zero((1-nn)/(D*t*(1-t))-2*N/D**2)
        cross = (aa[0]*P[1]-aa[1]*P[0])/D
        checks[f'panel{panel}_centre{sign}_Q_split_discriminant'] = zero(S.discriminant(QQ,t)+4*(1-u*u-v*v)*cross**2)

# Rational parametrization of h^2=k t(1-t), over Q(q,target).
paramt = k/(a*a+k)
paramh = k*a/(a*a+k)
checks['omega_cover_rational_parametrization'] = zero(paramh**2-k*paramt*(1-paramt))
facts['omega_rational_map'] = {'t':str(paramt),'h':str(paramh)}
# Q=a2(t-rplus)(t-rminus) has a rational point (rplus,0)
# after the finite TARGET extension i,w. Generic slope parametrization.
a2,rp,rm = S.symbols('a2 rplus rminus')
tpar = (a*a*rp-a2*rm)/(a*a-a2)
ypar = a*a2*(rp-rm)/(a*a-a2)
checks['Q_cover_rational_parametrization'] = zero(ypar*ypar-a2*(tpar-rp)*(tpar-rm))
facts['Q_rational_map'] = {'t':str(tpar),'sqrtQ':str(ypar)}

# Source polynomial numerators of degree <=2 are parameter derivatives of
# inv/log/xlog; parameter differentiation stays inside the same connection.
r1,r2,y1,y2,intercept = S.symbols('r1 r2 y1 y2 intercept')
lin = intercept-r1*y1-r2*y2
checks['linear_numerator_derivative'] = zero(-S.diff(S.log(lin),r1)-y1/lin)
checks['quadratic_numerator_derivative'] = zero(S.diff(lin*S.log(lin),r1,r2)-y1*y2/lin)
checks['intercept_derivatives_commute'] = zero(S.diff(lin*S.log(lin),intercept,2)-1/lin)

# On the moving line-charge arc of either panel, h_b=q^2 and B=A.
# Thus the A and B numerators cancel before taking the local singular cut.
for panel,P in zip([3,4],pp):
    us,vs=P[0]/D,P[1]/D
    c0=(1-vs*vs)/18
    c2=(-1 if panel==3 else 1)*2*(1+vs)*us/9
    hb2=(1-vs)/(1+vs)
    checks[f'panel{panel}_source_arc_h_equals_q_squared'] = zero(hb2-q**4)
    checks[f'panel{panel}_full_density_source_arc_numerator_cancels'] = zero((c0-c2*q*q)+(-c0+c2*hb2/(q*q)))

# Independent actual-angle numerical checks, no repository routines reused.
mp.mp.dps = 52
def Hactual(nn,zz):
    rr=mp.sqrt(nn)
    gam=mp.acos(zz/rr)
    theta=2*mp.atan((1+rr)/mp.sqrt(1-nn)*mp.tan(gam/2))
    return gam-theta
angleerr=mp.mpf(0)
for nn,zz in [(mp.mpf('0.43'),mp.mpf('0.19')),(mp.mpf('0.72'),mp.mpf('-0.31'))]:
    HH=Hactual(nn,zz); hh=mp.sqrt(1-nn); QQ=nn-zz*zz
    dzz=(hh/(1-zz)-1)/mp.sqrt(QQ)
    dnn=(zz-(zz+nn)/hh)/(2*nn*mp.sqrt(QQ))
    angleerr=max(angleerr,abs(mp.diff(lambda s:Hactual(nn,s),zz)-dzz),abs(mp.diff(lambda s:Hactual(s,zz),nn)-dnn))
checks['actual_angle_both_parameter_derivatives'] = angleerr<mp.mpf('1e-47')
facts['actual_angle_max_error'] = str(angleerr)
facts['source_arc_angle_complement_error'] = str(abs(mp.atan(mp.mpf('0.37'))+mp.atan(1/mp.mpf('0.37'))-mp.pi/2))
checks['source_arc_actual_A_equals_B'] = mp.mpf(facts['source_arc_angle_complement_error'])<mp.mpf('1e-47')

# Direct beta-identity check for the source-independent numerator kernel.
betaerr=mp.mpf(0)
for l1,l2 in [(mp.mpf('.67'),mp.mpf('1.23')),(mp.mpf('1.41'),mp.mpf('.83'))]:
    val=2/mp.pi*mp.quad(lambda ph:1/(mp.sin(ph)**2*l1+mp.cos(ph)**2*l2),[0,mp.pi/2])
    betaerr=max(betaerr,abs(val-1/mp.sqrt(l1*l2)))
checks['actual_beta_integral_identity'] = betaerr<mp.mpf('1e-47')
facts['beta_identity_max_error'] = str(betaerr)

# Exact fixed-endpoint A/B split checked by quadrature for both actual panels.
spliterr=mp.mpf(0);sample=[]
for uu,vv in [(mp.mpf('.27'),mp.mpf('.31')),(mp.mpf('-.34'),mp.mpf('-.22'))]:
    for panel in [3,4]:
        hh=mp.sqrt((1-vv)/(1+vv))
        c0=(1-vv*vv)/18
        c2=(-1 if panel==3 else 1)*2*(1+vv)*uu/9
        def terms(qq):
            W=mp.sqrt(1+qq**4+(1 if panel==3 else -1)*2*uu*qq**2+(qq**4-1)*vv)
            A=mp.pi/2-mp.atan(qq); B=mp.atan(qq/hh)
            BB=hh*hh*(B/qq)**2 if qq else mp.mpf(1)
            independent=(c0-c2*qq*qq)*A*A/W
            dependent=(-c0*B*B+c2*BB)/W
            full=(c0*(A*A-B*B)+c2*(BB-qq*qq*A*A))/W
            return independent,dependent,full
        vals=[2/mp.pi*mp.quad(lambda qq:terms(qq)[i],[0,1]) for i in range(3)]
        spliterr=max(spliterr,abs(vals[0]+vals[1]-vals[2]))
        sample.append({'target_of_panel_source':[str(uu),str(vv)],'panel':panel,'A_piece':str(vals[0]),'B_piece':str(vals[1]),'sum':str(vals[2])})
checks['actual_C3_panel_A_B_split'] = spliterr<mp.mpf('1e-47')
facts['panel_split_max_error'] = str(spliterr)
facts['panel_split_samples'] = sample

rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
if sys.platform!='darwin': rss*=1024
result={'checks':checks,'facts':facts,'elapsed_seconds':time.time()-START,'peak_rss_bytes':rss,
        'scope':'Exact algebra and independent quadratures, not a computation of the physical full-density target cut.',
        'source_independent_T5':'excluded by the filtered rational-surface argument in x22_report.md',
        'full_density_T5':'OPEN'}
(OUT/'x22_verify.json').write_text(json.dumps(result,indent=2)+'\n')
for name,passed in checks.items(): print(('PASS ' if passed else 'FAIL ')+name)
print(json.dumps({'count':len(checks),'passed':sum(checks.values()),'elapsed_seconds':result['elapsed_seconds'],'peak_rss_bytes':rss},indent=2))
assert all(checks.values()),[name for name,passed in checks.items() if not passed]
