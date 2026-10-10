"""PROOF (exact sympy for every algebraic step) of the elementary forms of the rational terms, x = arccos b in (0,pi),
rho = sqrt((1+b)/(1-b)) = cot(x/2),  L(t) = log((1+t)/(1-t)),  J(t) = log((1+b+(1-b)t^2)/(1+t^2)):
  R_rho = int_0^1 dt/(t^2+rho^2)          = (x/2) tan(x/2)
  L_rho = int_0^1 L t dt/(t^2+rho^2)      = x^2/4          (L1 = L_rho at x=pi/2 = pi^2/16, R1 = pi/4)
  L_reg = int_0^1 L 2t dt/((t^2+1)((1-b)t^2+1+b)) = (L1-L_rho)/b = (pi^2/4-x^2)/(4 cos x)
  JL    = int_0^1 J 2dt/(1-t^2)           = pi^2/8 - x^2/2
  J1    = int_0^1 J dt/(1+t^2)            = (pi/4) log((1+sin x)/2) + Ti2(tan(pi/4-x/2)),
  Ti2(tan u) = u log|tan u| + (Cl2(2u)+Cl2(pi-2u))/2   (|u|<pi/2).
Each: derivative in the parameter equals derivative of the closed form (exact), plus one value.  Numerical spot check at the end."""
import sympy as sp, mpmath as mp
t,x,rho,b,u=sp.symbols('t x rho b u',positive=True)
L=sp.log((1+t)/(1-t)); ok=True
def Z(e):
    e=sp.simplify(e)
    return e==0 or sp.simplify(sp.expand(e.rewrite(sp.exp)))==0
# --- L_rho: dF/drho = -2 rho int L t/(t^2+rho^2)^2 ; IBP with v(t) = -(1/2)(1/(t^2+rho^2) - 1/(1+rho^2)), v(1)=0
v=-sp.Rational(1,2)*(1/(t**2+rho**2)-1/(1+rho**2))
ok&=Z(sp.diff(v,t)-t/(t**2+rho**2)**2)                       # v' = t/(t^2+rho^2)^2
ok&=Z(sp.diff(L,t)*v + 1/((t**2+rho**2)*(1+rho**2)))          # -L' v = 1/((t^2+rho^2)(1+rho^2)) : rational integrand
# so int L t/(t^2+rho^2)^2 = (1/(1+rho^2)) int_0^1 dt/(t^2+rho^2) = atan(1/rho)/(rho(1+rho^2))  [boundary: L v -> 0 at t=1 (v=O(1-t)), L(0)=0]
dF=-2*rho*sp.atan(1/rho)/(rho*(1+rho**2))
# with rho = cot(x/2), atan(1/rho) = x/2 for x in (0,pi):
dFdx=(dF.subs(sp.atan(1/rho),x/2)*sp.diff(sp.cot(x/2),x)).subs(rho,sp.cot(x/2))
ok&=Z(dFdx-sp.diff(x**2/4,x)); print('L_rho derivative identity:',ok)
# value: L_rho -> 0 as rho -> oo (x -> 0): |L_rho| <= (1/rho^2) int_0^1 L t dt -> 0.
# --- R_rho
ok&=Z((sp.atan(1/rho)/rho).subs(sp.atan(1/rho),x/2).subs(rho,sp.cot(x/2))-x/2*sp.tan(x/2)); print('R_rho:',ok)
# --- L_reg partial fractions
ok&=Z(2*t/((t**2+1)*((1-b)*t**2+1+b))-(t/(t**2+1)-t/(t**2+(1+b)/(1-b)))/b); print('L_reg partial fractions:',ok)
# --- JL: d/db J = (1-t^2)/(1+b+(1-b)t^2), so d/db [J 2/(1-t^2)] = 2/((1-b)(t^2+rho^2)); int_0^1 = 2(x/2)/((1-b)rho) = x/sin x
J=sp.log((1+b+(1-b)*t**2)/(1+t**2))
ok&=Z(sp.diff(J,b)*2/(1-t**2)-2/((1-b)*(t**2+(1+b)/(1-b))))
ok&=Z(((1-b)*sp.sqrt((1+b)/(1-b))).subs(b,sp.cos(x))**2-sp.sin(x)**2)       # (1-b) rho = sin x (both positive)
# dJL/dx = (x/sin x) * db/dx = -x  = d/dx(pi^2/8 - x^2/2);  value 0 at b=0 (J == 0), RHS(pi/2) = 0
ok&=Z((x/sp.sin(x))*sp.diff(sp.cos(x),x)+x); print('JL:',ok)
# --- J1: d/db [J/(1+t^2)] = (1/b)[1/(1+t^2) - 1/((1-b)(t^2+rho^2))]
ok&=Z(sp.diff(J,b)/(1+t**2)-(1/(1+t**2)-1/((1-b)*(t**2+(1+b)/(1-b))))/b)
dJ1db=(sp.pi/4-x/(2*sp.sin(x)))/sp.cos(x)
dJ1dx=dJ1db*sp.diff(sp.cos(x),x)
Ti2=lambda w: sp.Function('Ti2')(w)
# derivative of RHS using Ti2'(w) = atan(w)/w, atan(tan(pi/4-x/2)) = pi/4-x/2 for x in (0,pi)
w=sp.tan(sp.pi/4-x/2)
dR=sp.diff(sp.pi/4*sp.log((1+sp.sin(x))/2),x)+(sp.pi/4-x/2)/w*sp.diff(w,x)
# the difference is linear in x and pi with trigonometric coefficients; rationalise with q = tan(x/4)
q=sp.Symbol('q',positive=True); Dd=sp.expand(dR-dJ1dx)
xs=sp.Symbol('xs'); Dd=Dd.subs(x,xs)  # placeholder
Dd=sp.expand(dR-dJ1dx)
X_,PI_=sp.symbols('X_ PI_')
def rat(e):  # replace free x by X_, then trig(x) by rational functions of q
    e=sp.expand_trig(sp.expand(e))
    e=e.subs({sp.sin(x):sp.sin(4*sp.atan(q)),sp.cos(x):sp.cos(4*sp.atan(q))})
    e=e.subs({sp.sin(x/2):2*q/(1+q**2),sp.cos(x/2):(1-q**2)/(1+q**2)})
    return e
Dd=Dd.subs(x,X_)  # free occurrences
Dd=Dd.subs({sp.sin(X_):sp.sin(x),sp.cos(X_):sp.cos(x),sp.tan(sp.pi/4-X_/2):sp.tan(sp.pi/4-x/2),sp.sin(X_/2):sp.sin(x/2),sp.cos(X_/2):sp.cos(x/2)})
Dd=sp.expand_trig(sp.expand(Dd))
Dd=sp.expand(Dd.subs({sp.sin(x):2*sp.tan(x/2)/(1+sp.tan(x/2)**2),sp.cos(x):(1-sp.tan(x/2)**2)/(1+sp.tan(x/2)**2)}).subs(sp.tan(x/2),q))
Dd=sp.cancel(sp.together(Dd.rewrite(sp.tan).subs(sp.tan(x/2),q)))
print('   J1 residual after rationalisation:',Dd)
ok&=(Dd==0); print('J1 derivative identity:',ok)
# value at x = pi/2 (b = 0): J == 0 so J1 = 0; RHS = (pi/4) log 1 + Ti2(0) = 0.
# --- Ti2(tan u) Clausen form: derivative  2u/sin(2u)  on both sides; Cl2'(th) = -log|2 sin(th/2)|
lhs=(u/sp.tan(u))*sp.diff(sp.tan(u),u)
# d/du RHS = [log tan u - log(2 sin u) + log(2 cos u)] + u tan'(u)/tan(u); the bracket is log(tan u cos u/sin u) = log 1 = 0
ok&=Z(sp.tan(u)*sp.cos(u)/sp.sin(u)-1)
ok&=Z(lhs-u*sp.diff(sp.tan(u),u)/sp.tan(u)); ok&=Z(lhs-2*u/sp.sin(2*u))
# value at u=0: Ti2(0)=0, Cl2(0)=Cl2(pi)=0.
print('ALL RATIONAL-TERM IDENTITIES PROVED' if ok else 'FAILURE')
# numerical spot check
mp.mp.dps=40
Cl2=lambda th: mp.clsin(2,th)
for xv in (mp.mpf('0.37'),mp.mpf('1.9'),mp.mpf('2.8')):
    bv=mp.cos(xv); rv=mp.cot(xv/2); Lf=lambda tt: mp.log((1+tt)/(1-tt)); Jf=lambda tt: mp.log((1+bv+(1-bv)*tt*tt)/(1+tt*tt))
    uu=mp.pi/4-xv/2
    e=[mp.quad(lambda tt: Lf(tt)*tt/(tt*tt+rv*rv),[0,1])-xv**2/4,
       mp.quad(lambda tt: Jf(tt)*2/(1-tt*tt),[0,1])-(mp.pi**2/8-xv**2/2),
       mp.quad(lambda tt: Jf(tt)/(1+tt*tt),[0,1])-(mp.pi/4*mp.log((1+mp.sin(xv))/2)+uu*mp.log(abs(mp.tan(uu)))+(Cl2(2*uu)+Cl2(mp.pi-2*uu))/2)]
    print('x=%s  numerical residuals:'%xv,[mp.nstr(abs(v),3) for v in e])
