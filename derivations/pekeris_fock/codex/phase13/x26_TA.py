import mpmath as m,pickle,json,resource,io
from pathlib import Path
m.mp.dps=50
p=Path(__file__).resolve().parent;r=p.parents[1]
class LegacyMpc:
 def __setstate__(self,state):self.state=state
class Compat(pickle.Unpickler):
 def find_class(self,module,name):
  if (module,name)==('mpmath.ctx_mp_python','mpc'):return LegacyMpc
  raise ValueError((module,name))
def realpart(t):
 sign,mant,exp,bc=t
 return (-1 if sign else 1)*m.mpf(int(mant,16))*m.power(2,exp)
def number(x):return m.mpc(realpart(x.state[0]),realpart(x.state[1]))
d=Compat(io.BytesIO((r/'claude10/TA.pkl').read_bytes())).load()
d={k:[[number(x) for x in row] for row in val] for k,val in d.items()}
A=m.matrix(d['T1'])-m.eye(5);B=m.matrix(d['T2'])-m.eye(5)
sv=lambda M:[m.nstr(x,10) for x in m.svd(M,compute_uv=False)]
res={'warning':'Saved 50-digit experiment, no rigorous error enclosure. Based line-A pair, not the X-strip pair.',
 'Nplus_square_norm':m.nstr(m.norm(A*A),10),'Nminus_square_norm':m.nstr(m.norm(B*B),10),
 'images_concatenated_singular_values':sv(m.matrix([[A[i,j] for j in range(5)]+[B[i,j] for j in range(5)] for i in range(5)])),
 'Nplus_Nminus_singular_values':sv(A*B),'Nminus_Nplus_Nminus_singular_values':sv(B*A*B),
 'alternating4_norm':m.nstr(m.norm((A*B)**2),10),'alternating6_norm':m.nstr(m.norm((A*B)**3),10),
 'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(p/'x26_TA.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))
