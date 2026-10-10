"""Symbolic H (Codex's closed auxiliary part + rr heads) on the same atoms, coefficients rational in T=tan(al/4), U=tan(ps/4)."""
import sympy as sp
from fractions import Fraction as Fr
from sym_classical import *
Ec, a21c = sp.symbols('E a21')
chih = -(3*pi+10-16*Gc)/(24*pi)
a = ca**2-sa**2; w = cp**2-sp_**2; d = r2*sp_; yv = -sp_; z = r2*cp; h = 2*ca*sa; sm = ca+sa
def logsmd(): return lam(1, 1, 1)+lam(1, -1, 5)-log2/2
def chi():
    s_, c_, g = h, a, 2*sp_*cp; sg = sm
    Lratio = 2*lam(-1, 1, 1)+2*lam(-1, -1, 5)-lam(-2, 2, 2)-lam(2, 2, 2)
    l5 = lam(2, -2, 2)+lam(2, -2, 6)+lam(2, 2, -2)-lam(2, 2, 2)-2*lam(-2, 2, 2)
    l6 = lam(2, -2, 2)+lam(2, 2, 2)-lam(2, -2, 6)-lam(2, 2, -2)
    ImD = cl(4, 4, -4)+cl(-4, -4, 12)-cl(-4, 4, 4)-cl(4, -4, 4)
    psi2 = ang(0, -1, Fr(1, 2))
    e = (g+s_-4*sg*d)/6+(-2*pi*w*logsmd()+pi*c_*Lratio+2*g*psi2+c_*psi2*l5+c_*AL*l6-c_*ImD)/(6*pi)+sg*d/3+chih*w
    return e
def bm():
    return (6*chih*w-6*chi()-2*w*logsmd()+2*ca*sa-2*sm*d+4*d*z*ang(0, -Fr(1, 2), Fr(1, 2))/pi)/a
def Rterms():
    A_, Y_, Hs, RC = sp.symbols('A_ Y_ Hs RC')
    LV, LC = sp.Function('LV')(A_, Y_), sp.Function('LC')(A_)
    hh = sp.sqrt(1-A_**2); t = (1-hh)/2; c = 2/(1+hh); v = Y_*sp.sqrt(c)
    m1 = -v*LV-1-v**3/2; m3 = -v**3*LV-sp.Rational(1, 3)-v**2; mm = -LV/v-v/2-v**3/4
    P_ = (m3-m1)/2; Q_ = (m3-2*m1+mm)/8
    R = (4*hh*P_+16*t*Q_)/sp.sqrt(c)+Y_*(2*t-LC)+Y_**3*(LC-c*t/2)
    dLV = {sp.Derivative(LV, A_): -sp.diff(v, A_)/(1-v), sp.Derivative(LV, Y_): -sp.diff(v, Y_)/(1-v)}
    dLC = {sp.Derivative(LC, A_): sp.diff(sp.log(c), A_)}
    Ra = sp.diff(R, A_).subs(dLV).subs(dLC); Ry = sp.diff(R, Y_).subs(dLV)
    expr = Z*r2/36*(A_*Ra/2-Y_*Ry-2*R)+Z**2/(6*r2)*(Y_*Ry+2*R)
    lv = lam(1, 1, 1)+lam(1, -1, 5)-lam(2, 0, 2); lcv = -2*lam(2, 0, 2)+2*log2
    expr = expr.subs({LV: lv, LC: lcv})
    expr = expr.subs(sp.sqrt(1-A_**2), Hs)
    expr = expr.subs(sp.sqrt(2)*sp.sqrt(1/(Hs+1)), RC).subs(sp.sqrt(2/(Hs+1)), RC).subs(sp.sqrt(1/(Hs+1)), RC/r2)
    assert not expr.has(A_) or True
    expr = expr.subs({Hs: h, RC: r2/sm, A_: a, Y_: yv})
    return expr
def S_(y):
    hh = (1-y*y)**sp.Rational(3, 2)
    return 4*sp.sqrt(pi)/3*hh+8/(3*sp.sqrt(pi))*(4*y**3/3-y+hh*ang(0, -Fr(1, 2), 0))
def W():
    M = r2*ca+r2*sa-2*cp*(1+2*ang(0, -Fr(1, 2), 0)/pi)
    asa = ang(-1, 0, Fr(1, 2)); bq = ang(-Fr(1, 2), 0, Fr(1, 4))
    Lr = 2*(2*bq*(lam(*form((-Fr(1, 2), 0, Fr(1, 4))))-lam(*form((-Fr(1, 2), 0, Fr(1, 4)), Fr(1, 2))))+cl(-4, 0, 4)+cl(4, 0, 4))
    Sy = 4*sp.sqrt(pi)/3*cp**3+8/(3*sp.sqrt(pi))*(4*yv**3/3-yv+cp**3*ang(0, -Fr(1, 2), 0))
    return Sy-sp.sqrt(pi)*a*a/(2*yv*yv)*(M+2*yv/pi*(1+h/a*asa))+sp.sqrt(pi)*a**3/(4*yv**3)*(bm()+Lr/pi)
def wRplus():
    u = (ca-sa)/r2; bb = 2*yv/a; v = -r2*sp_/sm
    LVs = lam(1, 1, 1)+lam(1, -1, 5)-lam(2, 0, 2)
    Fm = lambda m: -(LVs+sum(v**j/j for j in range(1, m+1)))/bb**(m+1)
    Lm = lambda m: u**(m+1)/(m+1)*LVs+bb*Fm(m+1)/(m+1)
    return 16/a*(Fm(0)-4*Fm(2)+4*Fm(4)-bb*(u*u/2-u**4/2)-bb**3*(u**4/4-u**6/2+u**8/4)-4/bb*(Lm(1)-2*Lm(3)))
def Hrr(which):
    if which == '+': v = -ca; lg = 2*lam(1, 0, 4)-log2; asv = ang(Fr(1, 2), 0, -Fr(1, 2)); sq = sa
    else: v = -sa; lg = 2*lam(1, 0, 2)-log2; asv = ang(-Fr(1, 2), 0, 0); sq = ca
    return sp.Rational(16, 3)-20*v*v+6*v**3+4*v*(3-5*v*v)*lg+w*(2*(5*v*v-3)/sq*(1+2*asv/pi)+12*v/pi)
def scalar_head():
    r1_, r2_ = ca, sa
    cshift = -sp.Rational(17, 72)+(24*Gc-31)/(36*pi); c = a21c+Z*cshift
    B = bm()
    aw = ang(0, -1, Fr(1, 2))/w; acd = ang(0, -Fr(1, 2), Fr(1, 2))
    angle_part = aw/12+w*aw**2/(12*pi)-d*z*acd/(3*pi)
    out = Ec*(Z/18*sm*(2+r1_*r2_)-(6-d*d)*d/72)-Z**3/18*sm*(1+5*r1_*r2_)
    out -= c*(Z/2*sm*w-(6-5*d*d)*d/12)
    out += (3-2*d*d)*d/72
    out += Z*(-(6-5*d*d)*d/36*logsmd()-a*d*B/24+5*d**3/(108*pi)-d**3/54-5*sm*d*d/36+(1+10*r1_*r2_)*d/72+sm**3/108+d*d*z*acd/(6*pi))
    out += Z**2*(sm*w*(logsmd()/6+1/(12*pi))+(r1_-r2_)*(1+sp.Rational(5, 2)*r1_*r2_)*B/9+sm*angle_part-29*d**3/216+sm*d*d/24+(17+20*r1_*r2_)*d/36-sm*(47-2*r1_*r2_)/216)
    return out
def H_all():
    return scalar_head()-5*Z*d**3/(54*pi)+Rterms()-Z*r2/(24*sp.sqrt(pi))*W()-Z**2/(192*r2)*wRplus()+Z**2/36*(Hrr('+')+Hrr('-'))
