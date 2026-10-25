"""Base-path/symmetry experiment, retaining its explicit identification caveat."""
from pathlib import Path
import json
p=Path(__file__).resolve().parent
s=(p/'x26_monodromy.py').read_text();exec(s[:s.index("A=load('op_e2_A.txt')")])
d=json.loads((p/'x26_monodromy.json').read_text())
M={k:m.matrix([[m.mpc(*x) for x in row] for row in a]) for k,a in d['matrices'].items()}
A=load('op_e2_A.txt');base=m.mpf(1)/5;l=m.sqrt(m.mpf(8)/9);sing=[l,-l,1,-1]
F=transport(A,[base,-base],sing);S=m.diag([1,-1,1,-1,1]);H=S*F
# H identifies the u-sign automorphism using its real nonsingular connector.
Uleft=H**-1*M['U_v_plus']*H;R=M['R_u_minus'];N=M['Q2']-m.eye(5)
res={'reflection_H_square':m.nstr(m.norm(H*H-m.eye(5)),8),'R_Uleft_commutator':m.nstr(m.norm(R*Uleft-Uleft*R),8),'cases':{}}
for tag,G in [('R_Uright',R*M['U_v_plus']),('R_Uright_inverse',R*M['U_v_plus']**-1),('R_Uleft',R*Uleft),('R_Uleft_inverse',R*(Uleft**-1)),('Uleft',Uleft)]:
 B=G*N*G**-1
 res['cases'][tag]={'image_singular_values':sv(m.matrix([[N[i,j] for j in range(5)]+[B[i,j] for j in range(5)] for i in range(5)])),'product':sv(N*B),'triple':sv(B*N*B)}
 print(tag,res['cases'][tag],flush=True)
res['H']=[[[m.nstr(m.re(H[i,j]),m.mp.dps),m.nstr(m.im(H[i,j]),m.mp.dps)] for j in range(5)] for i in range(5)]
(p/'x26_strip.json').write_text(json.dumps(res,indent=2))
