"""Exact S_e2 chart on the q=1 electron-panel boundary of the cap."""
import sys
sys.dont_write_bytecode=True
import resource,time,json
resource.setrlimit(resource.RLIMIT_CPU,(600,600))
from pathlib import Path
import sympy as S
START=time.time();OUT=Path(__file__).resolve().parent
xi,tau,T,u,v=S.symbols('xi tau T u v');xx=S.Matrix([u,v]);a=S.Matrix([0,1]);D=1+xi**4;N=(1+xi*xi)**2
c=(1-xi*xi)/(1+xi*xi);tau_minus=-2*xi*xi/D
maptau=tau_minus+(1-tau_minus)*T;jac=S.factor(1-tau_minus)
checks={};facts={'panels':[]}
def zero(x):return S.cancel(x)==0
def test(name,b):
 checks[name]=bool(b);print(('PASS ' if b else 'FAIL ')+name,flush=True)
 if not b:raise AssertionError(name)
for panel,sign in [(1,-1),(2,1)]:
 d=S.Matrix([sign,0]);R=(1-tau)*c*d+tau*a
 pminus=S.simplify(R.subs(tau,tau_minus))
 pp=S.Matrix([sign*(1-xi**4),-2*xi*xi])/D
 test(f'panel{panel}_second_unit_intersection',all(zero(x) for x in pminus-pp))
 test(f'panel{panel}_unit_pminus',zero(pp.dot(pp)-1))
 Rnew=T*a+(1-T)*pp
 test(f'panel{panel}_affine_chord_map',all(zero(x) for x in R.subs(tau,maptau)-Rnew))
 nn=S.factor(Rnew.dot(Rnew));zz=Rnew.dot(xx);Q0=S.factor(nn-zz*zz)
 test(f'panel{panel}_H_beta',zero(1-nn-2*N*T*(1-T)/D))
 # Original q=1 contributes 1/sqrt(2). The old t beta integral is pi.
 # Current radical two-form / omega_K3 has ratio i*M on the chosen Q0 root.
 M=xi*N/D*T*(1-T)
 ratio_squared=S.cancel(xi*xi/((1+xi*xi)**2*2)*(1-nn)/(-Q0)*jac*jac*D*T*(1-T)*Q0)
 test(f'panel{panel}_form_ratio_square',zero(ratio_squared+M*M))
 test(f'panel{panel}_odd_multiplier',zero(M.subs(xi,-xi)+M))
 B=S.factor(a.dot(D*pp)-(a.dot(xx))*(D*pp).dot(xx))
 A=1-v*v;C=S.factor(D*D-((D*pp).dot(xx))**2)
 facts['panels'].append({'panel':panel,'density_direction':[str(x) for x in d],
  'pminus':[str(S.factor(x)) for x in pp],
  'R':[str(S.factor(x)) for x in Rnew], 'n':str(nn),'z':str(zz),'Q0':str(Q0),
  'weierstrass':'Y^2=X^3+2*B*X^2+A*C*X','A':str(A),'B':str(B),'C':str(C),
  'standard_identification':'standard q=i*xi' if panel==1 else 'P=-P_standard(xi), then X->-X,Y->iY',
  'physical_chain':'0<xi<1; 2*xi^2/(1+xi^2)^2<T<1; retain z>1 for the target cap',
  'face_density':'log((1+xi)/(1-xi))*C_actual(n,z,R,u,v,q=1)*i*M*omega_K3, plus the rational cut character',
  'M':str(M)})
# The rational quotient under xi->-xi has no transcendental H2 (x19/x20).
# A log-weighted integral must NOT be projected merely by this parity.
facts.update({'variables':['xi','tau','T','u','v'],'tau_minus':str(tau_minus),
 'map':'tau=tau_minus+(1-tau_minus)*T','jacobian':str(jac),
 'T_min':str(S.factor(-tau_minus/(1-tau_minus))),
 'omega_K3':'d_xi wedge d_T/sqrt((1+xi^4)*T*(1-T)*(n-z^2))',
 'measure_identity':'[xi/(sqrt(2)*(1+xi^2))]*sqrt((1-n)/(z^2-n))*d_xi wedge d_tau = i*M*omega_K3; choose sqrt(n-z^2)=i*sqrt(z^2-n) on the cap',
 'old_beta_integral':'Integral(dt/sqrt(t*(1-t)),0,1)=pi; R is independent of old t on this face',
 'known_input':'This is the same bare S_e2 K3 family, after the explicitly stated complex base changes. Rank five/irreducibility are CITED x20 inputs, not recomputed here.',
 'parity_warning':'M is odd in xi and omega_K3 is odd under xi->-xi. Without the density log this component descends to the rational quotient. The log is odd too and its derivative can couple the graded pieces. Parity does not exclude a middle T5 factor of the log extension.',
 'physical_warning':'q=1 is a boundary of the four-dimensional cap, not by itself a physical target variation or closed Betti cycle. No boundary/Gysin coefficient from the full relative mapping cone and no nonzero T5 pairing is asserted.',
 'next_obstruction':'Need a characteristic-zero certificate for the log-extension transgression and relative q=1/cap boundary pairing, including the other faces and electron sum. Identifying this boundary surface alone does not supply that map.'})
test('physical_T_lower_endpoint',zero((-tau_minus/(1-tau_minus))-2*xi*xi/N))
facts['checks']=checks;facts['elapsed_seconds']=time.time()-START
rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss;facts['peak_rss_bytes']=rss if sys.platform=='darwin' else rss*1024
(OUT/'x24_q1.json').write_text(json.dumps(facts,indent=2)+'\n')
(OUT/'x24_q1.txt').write_text('# Exact q=1 boundary chart to S_e2. Physical framing remains OPEN.\n'+json.dumps(facts,indent=2)+'\n')
print(json.dumps({'checks':len(checks),'elapsed_seconds':facts['elapsed_seconds'],'peak_rss_bytes':facts['peak_rss_bytes']},indent=2))
