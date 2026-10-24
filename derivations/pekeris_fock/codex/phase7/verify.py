"""Small serial checks of exported matrices, physical pullbacks and interfaces."""
from pathlib import Path
import json,sys,hashlib
import sympy as S
import mpmath as m
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
sys.path.insert(0,str(OUT.parent/'phase3'))
import panel_connection as P
from classical_reference import L,T,lc
m.mp.dps=45

def vector(names):return [m.mpf((i*7)%23-11)/19 for i in range(len(names))]
def eval_action(matrix,variables,values,X):
 symbols=S.symbols(' '.join(variables));loc=dict(zip(variables,symbols));loc['I']=S.I
 exprs=[S.sympify(row[2],locals=loc) for row in matrix['entries']]
 f=S.lambdify(symbols,exprs,'mpmath',cse=True)
 out=[m.mpc(0) for _ in X]
 for (i,j,_,_),c in zip(matrix['entries'],f(*values)):out[i]+=c*X[j]
 return out

def padded(names,X,family):
 lookup=dict(zip(names,X));js=['00','10','20'] if family=='A' else ['00','01','02'];out=[]
 for fam in ['J','K']+['J'+str(k) for k in range(8)]+['FP','E','F0']:
  out.append([lookup.get(f'{fam}_{j}',m.mpf(0) if fam=='F0' and j!='00' else m.mpf(2)/7) for j in js])
 return out

def run():
 d=json.loads((OUT/'connections.json').read_text());names=d['atoms'];X=vector(names);worst=0;checks=0
 for b,Y in [(m.mpf('.3'),m.mpf('-.6')),(m.mpf('-.27'),m.mpf('-.7'))]:
  side='plus' if b>0 else 'minus';h,zs,u=P.geometry(b,Y);z=zs[0];expected=[]
  for axis,axisname in [(0,'b'),(1,'y')]:
   action=eval_action(d['panels'][side+'_'+axisname],['b','y','h','z','u'],[b,Y,h,z,u],X)
   ref={}
   for family in ['A','D']:
    js=['00','10','20'] if family=='A' else ['00','01','02']
    rows=P.rhs(b,Y,axis,family,padded(names,X,family))
    for fam,row in zip(['J','K']+['J'+str(k) for k in range(8)]+['FP','E','F0'],rows):
     for jet,value in zip(js,row):
      name=f'{fam}_{jet}'
      if name in names:ref[name]=value
   error=max(abs(a-ref[name]) for a,name in zip(action,names));worst=max(worst,error);checks+=55
   assert error<m.mpf('1e-37'),(side,axis,error)
   expected.append(action)
  F=lambda H,Z:-(Z*Z-1)*(Z*Z-H*H)/(2*Z*Z*(1+H*H))
  bh=4*h/(1+h*h)**2;Yh=m.diff(lambda H:F(H,z),h)/(2*Y);Yz=m.diff(lambda Z:F(h,Z),z)/(2*Y)
  for ax,coefs in [('h',[bh,Yh]),('z',[0,Yz])]:
   action=eval_action(d['cover'][side+'_'+ax],['h','z'],[h,z],X)
   target=[coefs[0]*a+coefs[1]*bb for a,bb in zip(*expected)]
   error=max(abs(a-bb) for a,bb in zip(action,target));worst=max(worst,error);checks+=55
   assert error<m.mpf('1e-37'),(side,ax,error)
 print('PASS panel/cover actions and pullbacks:',checks,'components; max error',m.nstr(worst,5),flush=True)
 def vals(a,p):
  out=[]
  for sign in [1,-1]:
   w=(a+sign*p)/4;out.extend(L(w+mm*m.pi/4) for mm in [-1,0,1]);out.append(T(w))
  out.extend([L(p/2),L((a+m.pi/2)/2),L((m.pi/2-a)/2),lc((a+p)/4),lc((a-p)/4),a,p,a*a,p*p,m.mpf(1)])
  return out
 a=m.mpf('-.4');p=m.mpf('-1.3');t=m.exp(1j*a/2);q=m.exp(1j*p/2);X=vals(a,p);err=0
 for axis,an,var in [(0,'dt',t),(1,'dq',q)]:
  got=eval_action(d['classical'][an],['t','q'],[t,q],X)
  for k,v in enumerate(got):
   expected=m.diff(lambda arg:vals(arg,p)[k] if axis==0 else vals(a,arg)[k],a if axis==0 else p)
   err=max(err,abs(v*(1j*var/2)-expected))
 assert err<m.mpf('1e-38');print('PASS 36 actual classical derivatives; max error',m.nstr(err,5),flush=True)
 Y=m.mpf('-.6')
 def integral(b):
  h=m.sqrt((1+b)/(1-b));u=min(m.mpf(1),h)
  return m.quad(lambda s:1/m.sqrt(1+b+(4*Y*Y-2)*s*s+(1-b)*s**4),[0,u])
 step=m.mpf('1e-12');mid=integral(0);jump=(integral(step)-2*mid+integral(-step))/step
 predicted=-1/(2*abs(Y));assert abs(jump-predicted)<m.mpf('1e-10')
 print('PASS actual endpoint derivative jump:',m.nstr(jump,20),'predicted',m.nstr(predicted,20),flush=True)
 rec=json.loads((OUT/'angle_recurrence.json').read_text())
 assert rec['stats']=={'diag':{'count':563,'shell_range':[-7,7]},'lower1':{'count':415,'shell_range':[-7,7]},'lower2':{'count':251,'shell_range':[-6,6]}}
 print('PASS saved independent angle-subblock counts',flush=True)
 inputs=['codex/phase3/panel_connection.py','codex/phase3/classical_reference.py','codex/codex_derive.py','codex/phase5/goal2_combined_exact.py','codex/phase5/goal2_pieces_exact.json']
 (OUT/'input_hashes.json').write_text(json.dumps({p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},indent=2)+'\n')
 print('PASS all serial verification checks',flush=True)
if __name__=='__main__':run()
