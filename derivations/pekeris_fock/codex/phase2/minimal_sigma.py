#!/usr/bin/env python3
"""Two replacement blocks for the minimal five-atom log(sigma) system.

Other one-target blocks are pi times the corresponding six-atom block;
other rows are unchanged, Slog row/column is deleted, L is renamed Ltilde.
"""
import sys
sys.dont_write_bytecode=True
import ast,json
from pathlib import Path
import sympy as S
import order2 as O
import recurrence as R
OUT=Path(__file__).resolve().parent

def main():
    raw=json.loads((OUT/'cleared_operators.json').read_text())
    result={};operators={};stats={}
    for kind in ('S2P','C1'):
        left={ast.literal_eval(d):R.parse(f) for d,f in raw[kind]['(one,L)'].items()}
        right={ast.literal_eval(d):R.parse(f) for d,f in raw[kind]['(one,Slog)'].items()}
        op={d:S.expand(S.pi*left.get(d,0)-(S.pi-2)*right.get(d,0)) for d in set(left)|set(right)}
        op={d:f for d,f in op.items() if f!=0}
        stencil=R.rec(op)
        fixed={d:f.compose(R.energy,-R.ee**2) for d,f in stencil.items()}
        fixed={d:f for d,f in fixed.items() if f}
        result[kind]={'(one,Ltilde)':{str(d):str(f).replace('lam','lambda').replace('PI','pi') for d,f in sorted(stencil.items())}}
        operators[kind]={'(one,Ltilde)':{str(d):str(f) for d,f in op.items()}}
        stats[kind]={'(one,Ltilde)':{'general_count':len(stencil),'fixed_count':len(fixed),
            'general_shell_range':[min(map(sum,stencil)),max(map(sum,stencil))],
            'fixed_shell_range':[min(map(sum,fixed)),max(map(sum,fixed))]}}
        print(kind,stats[kind],flush=True)
    (OUT/'minimal_sigma_corrections.json').write_text(json.dumps(result,indent=2)+'\n')
    (OUT/'minimal_sigma_operators.json').write_text(json.dumps(operators,indent=2)+'\n')
    (OUT/'minimal_sigma_counts.json').write_text(json.dumps(stats,indent=2)+'\n')

if __name__=='__main__':main()
