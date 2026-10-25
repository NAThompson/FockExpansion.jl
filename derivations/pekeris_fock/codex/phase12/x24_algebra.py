"""Exact algebra for the physical B-density pilot; own code, no imported inputs."""
import sys
sys.dont_write_bytecode = True
import resource
resource.setrlimit(resource.RLIMIT_CPU, (600, 600))
import json, time
from pathlib import Path
import sympy as S
START=time.time(); OUT=Path(__file__).resolve().parent
n,z,b,q,xi,t,tau,eta,r,s,u,v,Y1,Y2=S.symbols('n z b q xi t tau eta r s u v Y1 Y2')
Q=n-z*z; om=1-n; ell=1-z
checks={}; facts={}
def zero(x):return S.cancel(x)==0
def test(name, cond):
    checks[name]=bool(cond)
    print(('PASS ' if cond else 'FAIL ')+name,flush=True)
    if not cond:raise AssertionError(name)
def strings(mat):return [[str(S.factor(x)) for x in row] for row in mat.tolist()]

# Density connection: column state F=(B^2,hB,1); h^2=(1+b)/(1-b).
Lq=1+b+(1-b)*q*q
Db=S.Matrix([[0,-2*q/((1+b)*Lq),0],[0,1/(1-b*b),-q/((1-b)*Lq)],[0,0,0]])
Dq=S.Matrix([[0,2*(1-b)/Lq,0],[0,0,(1+b)/Lq],[0,0,0]])
test('density_rank3_flat',all(zero(x) for x in Db.diff(q)-Dq.diff(b)+Db*Dq-Dq*Db))
facts['density']={'variables':['b','q'],'basis':['B^2','h*B','1'],'definitions':{'h^2':'(1+b)/(1-b)','B':'atan(q/h)'},'d_b':strings(Db),'d_q':strings(Dq),'graded_characters':['1','h','1']}

# Exact reduction of the two density parameters to one. The elementary
# integration-by-parts proof uses the derivative below, not an illegal PV split.
a=S.symbols('a')
log_arg=(1+a*a*xi*xi)/(1+a*a)
integrand=-S.diff(S.log(log_arg),a)/(1-xi*xi)
test('one_density_parameter_derivative',zero(integrand-2*a/((1+a*a)*(1+a*a*xi*xi))))
Lr=1+b+(1-b)*r*r; Ls=1+b+(1-b)*s*s
partial_num=(1-r*r)/Lr-(1-s*s)/Ls
test('density_partial_fraction_identity',zero(partial_num-2*(s*s-r*r)/(Lr*Ls)))
# The subtracted numerator at r=0 is essential for endpoint/Fubini validity.
test('density_subtracted_endpoint_identity',zero((1-r*r)/Lr-1/(1+b)+2*r*r/((1+b)*Lr)))
facts['density_reduction']={'exact_identity':'atan(q/h)^2 = (1-b)*q^2*Integral(xi*log((1+xi)/(1-xi))/(1+b+(1-b)*q^2*xi^2), (xi,0,1))',
 'proof':'Integrate by parts after subtracting log(1+a^2), a=q/h. Differentiate in a: 2*atan(a)/(1+a^2); value at a=0 is zero.',
 'endpoint_warning':'Do not split the two-term partial fraction into unsubtracted principal values: the corner r=s=0 carries a missing contribution.'}

# Zonal column Z=(K=h*J,log(1-z),1,h), h^2=1-n.
Az=S.Matrix([[z/Q,0,om/(ell*Q),-1/Q],[0,0,-1/ell,0],[0,0,0,0],[0,0,0,0]])
An=S.Matrix([[-1/(2*om)-1/(2*Q),0,-(z+n)/(2*n*Q),z/(2*n*Q)],
             [0,0,0,0],[0,0,0,0],[0,0,0,-1/(2*om)]])
test('zonal_rank4_flat',all(zero(x) for x in Az.diff(n)-An.diff(z)+Az*An-An*Az))
facts['zonal']={'basis':['K=sqrt(1-n)*J','log(1-z)','1','sqrt(1-n)'],
 'J':'(gamma-theta)/sqrt(n-z^2)', 'gamma':'acos(z/sqrt(n))',
 'theta':'acos((z-n)/(sqrt(n)*(1-z))) with physical continuous branches',
 'd_n':strings(An),'d_z':strings(Az),
 'graded_characters':['sqrt((1-n)/(n-z^2))','1','1','sqrt(1-n)'],
 'qualification':'All inverse rows below are modulo target polynomials of degree <=2; the actual regular pure inverse fixes these. Polynomials have no target cut.'}
R3=z*(4*z*z-3*n)
C1=-(n-4)*R3/(12*n**3)
A1=(4*z*z-n)*(3*n-4)/(12*n**3)
Cl=om*R3/(12*n**3)
Al=-(3*n*n*z*z+3*n*n-12*n*z*z-2*n+8*z*z)/(24*n**3)
Cx=om*om*R3/(60*n**3)
Ax=-(-6*n**3*z+15*n*n*z*z+5*n*n-20*n*z*z-2*n+8*z*z)/(120*n**3)
U1=S.Matrix([[C1,A1,0,0]]); Ul=S.Matrix([[Cl,Al,0,0]]); Ux=S.Matrix([[Cx,Ax,0,0]])
def cov(row,var):return (row.diff(var)+row*(An if var==n else Az)).applyfunc(S.cancel)
def euler(row):return (2*n*cov(row,n)+z*cov(row,z)).applyfunc(S.cancel)
U2=(U1+euler(U1)).applyfunc(S.cancel)
U3=(U2+euler(U2)/2).applyfunc(S.cancel)
rows={'inverse1':U1,'inverse2':U2,'inverse3':U3,'log':Ul,'ell_log':Ux}
facts['inverse_rows']={k:strings(m)[0] for k,m in rows.items()}
# Source identities ensure that the finite numerator prescription is exact.
x1,x2,rr1,rr2=S.symbols('x1 x2 R1 R2'); ll=1-rr1*x1-rr2*x2
for degree,base,factor in [(0,ll**-2,1),(1,ll**-1,1),(2,S.log(ll),-1),(3,ll*S.log(ll),1)]:
    val=S.diff(base,rr1,degree)*factor
    test('four_parameter_numerator_order_'+str(degree),zero(val-x1**degree/ll**2))
for degree,base,factor in [(0,ll**-3,1),(1,ll**-2,S.Rational(1,2)),(2,ll**-1,S.Rational(1,2)),(3,S.log(ll),-S.Rational(1,2)),(4,ll*S.log(ll),S.Rational(1,2))]:
    test('six_parameter_numerator_order_'+str(degree),zero(factor*S.diff(base,rr1,degree)-x1**degree/ll**3))

# Actual rational charge coefficients; sqrt(2) is factored out of panels 1,2.
D=1+q**4; yy=S.Matrix([Y1,Y2]); xx=S.Matrix([u,v])
P={1:S.Matrix([q**4-1,2*q*q]),2:S.Matrix([1-q**4,2*q*q]),
   3:S.Matrix([-2*q*q,1-q**4]),4:S.Matrix([2*q*q,1-q**4])}
panels=[]
for j in range(1,5):
    bj=Y1 if j==1 else -Y1 if j==2 else -Y2
    d=S.Matrix([-1,0]) if j==1 else S.Matrix([1,0]) if j==2 else S.Matrix([0,1])
    test('panel'+str(j)+'_b_direction',zero(bj+d.dot(yy)))
    test('panel'+str(j)+'_unit_arc',zero(P[j].dot(P[j])-D*D))
    charges=[1,2] if j in (1,2) else [2]
    for charge in charges:
        if j in (1,2):
            c0=(bj*bj-1)/72 if charge==1 else (2-bj*bj)/64
            c2=(1-bj)*(9*bj-8*Y2)/144 if charge==1 else 5*(1-bj)*Y2/32
            sector=charge; centres=[S.Matrix([0,1])]; sign=-1
        else:
            c0=(1-Y2*Y2)/18
            c2=(-1 if j==3 else 1)*2*(1+Y2)*Y1/9
            sector=3; centres=[S.Matrix([1,0]),S.Matrix([-1,0])];sign=1
        F4=S.factor(-q*q*c0*(1-bj)+c2*(1+bj))
        F6=S.factor(-c0*(1-bj*bj)+c2*(1+bj)**2/q**2)
        test(f'panel{j}_charge{charge}_density_numerator',zero((1-bj)*(-c0+c2*(1+bj)/((1-bj)*q*q))*q*q-F4))
        test(f'panel{j}_charge{charge}_four_numerator_degree',S.Poly(F4,Y1,Y2).total_degree()<=3)
        test(f'panel{j}_charge{charge}_six_numerator_degree',S.Poly(F6,Y1,Y2).total_degree()<=4)
        for aa in centres:
            panels.append({'panel':j,'psi30_charge':charge,'psi40_sector':sector,'a':[str(x) for x in aa],
                           'd':[str(x) for x in d],'P':[str(x) for x in P[j]],'b':str(bj),
                           'c0_rational':str(c0),'c2_rational':str(c2),
                           'four_numerator':str(F4),'six_numerator':str(F6),
                           'four_prefactor':str(sign*4)+'*sqrt(2)/pi^2',
                           'six_prefactor':str(sign*8)+'*sqrt(2)/pi^2'})
facts['physical_terms']=panels
# Electrons: u->-u sends panels 1<->2,3<->4 and nuclear centres reverse.
for charge in (1,2):
    f=[p for p in panels if p['panel']==1 and p['psi30_charge']==charge][0]
    g=[p for p in panels if p['panel']==2 and p['psi30_charge']==charge][0]
    test('electron_exchange_electron_sector_'+str(charge),zero(S.sympify(f['four_numerator'],locals=globals()).subs(Y1,-Y1)-S.sympify(g['four_numerator'],locals=globals())))
f=[p for p in panels if p['panel']==3][0];g=[p for p in panels if p['panel']==4][0]
test('electron_exchange_nuclear_sector',zero(S.sympify(f['four_numerator'],locals=globals()).subs(Y1,-Y1)-S.sympify(g['four_numerator'],locals=globals())))

# Simplex Jacobians, measured relative to d lambda_a,d lambda_p[,d lambda_r].
la=tau*t;lp=tau*(1-t);lr=(1-tau)*eta
jac6=S.det(S.Matrix([la,lp,lr]).jacobian([t,tau,eta]))
jac4=S.det(S.Matrix([la,lp]).jacobian([t,tau]))
test('six_simplex_jacobian',zero(jac6-tau*(1-tau)))
test('four_simplex_jacobian',zero(jac4-tau))
facts['chain_six']={'domain':'0<q,xi_r,xi_s,t,tau,eta<1; r=q*xi_r,s=q*xi_s',
 'R':'tau*(t*a+(1-t)*P/D)+(1-tau)*(eta*kappa(r)+(1-eta)*kappa(s))*d',
 'kappa':'kappa(s)=(1-s^2)/(1+s^2)',
 'measure':'q^2*(1-tau)/(sqrt(D*t*(1-t))*(1+r^2)*(1+s^2))',
 'source_kernel':'F6(y)/(1-R.y)^3',
 'faces':['q=0','q=1','xi_r=0','xi_r=1','xi_s=0','xi_s=1','t=0','t=1','tau=0','tau=1','eta=0','eta=1']}
facts['chain_four']={'domain':'0<q,xi,t,tau<1',
 'R':'tau*(t*a+(1-t)*P/D)+(1-tau)*kappa(q*xi)*d',
 'measure':'xi*log((1+xi)/(1-xi))/(sqrt(D*t*(1-t))*(1+q^2*xi^2))',
 'source_kernel':'F4(y)/(1-R.y)^2',
 'faces':['q=0','q=1','xi=0','xi=1','t=0','t=1','tau=0','tau=1']}

# Exact univariate chain-cut incidence: z=z_d+tau*(z_0-z_d).
z0,zd,nn0=S.symbols('z0 zd n0')
tauc=(1-zd)/(z0-zd)
test('cap_simplex_intersection',zero(1-(zd+tauc*(z0-zd))))
facts['candidate_target_cut']={'path':'base (u,v)=(1/5,1/3); vary v to 6/5 in upper/lower half planes, fix u=1/5. This is a prescription, not certified global chain transport.',
 'slice':'z_d=kappa(q*xi)*(d.x), z_0=(t*a+(1-t)*P/D).x',
 'cap_boundary':'tau_c=(1-z_d)/(z_0-z_d), when 0<tau_c<1',
 'delta_face_factor':'delta(z-1) d_tau -> 1/abs(z_0-z_d), if transversal',
 'status':'Local jump and cap incidence verified; endpoint regularization, complete global winding chain and T5 pushforward remain OPEN.'}

# Cut of the regular zonal solutions. Let sroot=sqrt(z^2-n)>0, z>1.
# Upper minus lower continuation in z is -2*pi*i*(A+C*h/sroot).
hh,ss=S.symbols('hh ss')
def reduce_roots(expr):
    num,den=S.fraction(S.cancel(expr))
    # Rationalize each square-root variable before reduction.
    expr=S.cancel(expr.subs(ss**2,z*z-n).subs(hh**2,om))
    num,den=S.fraction(expr)
    num=S.rem(S.Poly(num,ss),S.Poly(ss*ss-z*z+n,ss)).as_expr()
    den=S.rem(S.Poly(den,ss),S.Poly(ss*ss-z*z+n,ss)).as_expr()
    # comparison after denominator-clearing only
    return S.factor(num),S.factor(den)
def derivative_rad(expr,var):
    return S.diff(expr,var)+S.diff(expr,hh)*(-1/(2*hh) if var==n else 0)+S.diff(expr,ss)*(-1/(2*ss) if var==n else z/ss)
def E_rad(expr):return 2*n*derivative_rad(expr,n)+z*derivative_rad(expr,z)
G=A1+C1*hh/ss
G2=G+E_rad(G); G3=G2+E_rad(G2)/2
sub1={z:S.Integer(1),ss:hh}
g_at=S.factor(G.subs(sub1).subs(hh**2,om))
test('inverse1_jump_zero_at_pole',zero(g_at))
gz_at=S.factor(derivative_rad(G,z).subs(sub1).subs(hh**2,om))
test('inverse1_jump_normal_derivative',zero(gz_at+1/(4*om)))
facts['local_jump']={'orientation':'upper z+i0 minus lower z-i0; 0<n<1; z>1; sqrt(z^2-n)>0',
 'U1_jump':'-2*pi*i*Theta(z-1)*G1',
 'G1':'A1(n,z)+C1(n,z)*sqrt(1-n)/sqrt(z^2-n)',
 'U2_jump':'-2*pi*i*Theta(z-1)*(1+E)*G1',
 'U3_jump':'-2*pi*i*(Theta(z-1)*(1+E/2)*(1+E)*G1-delta(z-1)/(8*(1-n)))',
 'E':'2*n*d_n+z*d_z',
 'log_jump':'-2*pi*i*Theta(z-1)*(Alog+Clog*sqrt(1-n)/sqrt(z^2-n))',
 'ell_log_jump':'-2*pi*i*Theta(z-1)*(Axlog+Cxlog*sqrt(1-n)/sqrt(z^2-n))',
 'polynomial_numerators':'Apply the exported parameter derivative recipes DISTRIBUTIONALLY. Moving-cap faces and delta terms must not be dropped.',
 'scope':'These are local kernel jumps at fixed parameters, not a certified global physical cut or a five-state T5 vector.'}
# Homogeneous jump checks use ss^2=z^2-n, hh^2=1-n.
for label,AA,CC in [('inverse1',A1,C1),('log',Al,Cl),('ell_log',Ax,Cx)]:
    GG=AA+CC*hh/ss
    dd=derivative_rad(derivative_rad(GG,z),z)
    residual=S.cancel((n-z*z)*dd-3*z*derivative_rad(GG,z)+8*GG)
    nn,dden=S.fraction(residual)
    reduced=S.rem(S.Poly(nn,ss),S.Poly(ss*ss-z*z+n,ss)).as_expr()
    # Only inverse1 has single-valued source, so its jump is homogeneous.
    if label=='inverse1':test('inverse1_local_jump_homogeneous',zero(reduced))

# Defining equations exported unexpanded for finite-field verification.
facts['geometry']={'variables':['q','xi','t','tau','u','v'],
 'D':'1+q^4','p':'P/D','c':'(1-q^2*xi^2)/(1+q^2*xi^2)',
 'R':'tau*(t*a+(1-t)*p)+(1-tau)*c*d','n':'R1^2+R2^2','z':'R1*u+R2*v',
 'top_cover':'W^2=D*t*(1-t)*(1-n)*(n-z^2)',
 'top_character':'mu*sqrt((1-n)/(n-z^2)) with mu=1/sqrt(D*t*(1-t)); rational factors omitted',
 'warning':'The cover is FOUR-dimensional. It has not been identified with S_e2 or reduced to its T5 quotient.'}
# Generic endpoint radicals for every panel/centre; avoid expanding 6 variables.
for idx,item in enumerate(panels):
    aa=S.Matrix([S.sympify(x) for x in item['a']]);dd=S.Matrix([S.sympify(x) for x in item['d']]);pp=P[item['panel']]/D
    N=S.factor(D-aa.dot(P[item['panel']]))
    rr0=t*aa+(1-t)*pp;nn=S.factor(rr0.dot(rr0))
    test('term'+str(idx)+'_tau1_radical_cancellation',zero((1-nn)/(D*t*(1-t))-2*N/D**2))
    cc=(1-q*q*xi*xi)/(1+q*q*xi*xi)
    test('term'+str(idx)+'_tau0_h_square',zero(1-cc*cc-(2*q*xi/(1+q*q*xi*xi))**2))

facts['checks']=checks
facts['elapsed_seconds']=time.time()-START
rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
facts['peak_rss_bytes']=rss if sys.platform=='darwin' else rss*1024
(OUT/'x24_exact.json').write_text(json.dumps(facts,indent=2)+'\n')
lines=['# x24 exact rational handoff, generated by x24_algebra.py', '# All rows act on COLUMN bases; dF=A*F. Rational expressions use ** in JSON.',json.dumps(facts,indent=2)]
(OUT/'x24_exact.txt').write_text('\n'.join(lines)+'\n')
print(json.dumps({'passed':len(checks),'elapsed_seconds':facts['elapsed_seconds'],'peak_rss_bytes':facts['peak_rss_bytes']},indent=2))
