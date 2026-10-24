#!/usr/bin/env python3
"""Independent exact action checks for saved expanded stencils.

Run after recurrence.py. No regeneration or use of its finite1 kernel.
"""
import sys
sys.dont_write_bytecode=True
sys.setrecursionlimit(20000)  # expanded polynomial sums have deep additive ASTs
import ast,json,math,random
from fractions import Fraction as Q
from pathlib import Path
import sympy as S
import order2 as O
import recurrence as R
OUT=Path(__file__).resolve().parent

class ExactDivision(ast.NodeTransformer):
    def visit_BinOp(self,node):
        self.generic_visit(node)
        if isinstance(node.op,ast.Div):
            return ast.copy_location(ast.Call(func=ast.Name(id='Q',ctx=ast.Load()),
                       args=[node.left,node.right],keywords=[]),node)
        return node

def compile_coefficient(text):
    tree=ast.parse(text.replace('lambda','lam'),mode='eval')
    permitted=(ast.Expression,ast.BinOp,ast.UnaryOp,ast.Add,ast.Sub,ast.Mult,
               ast.Div,ast.Pow,ast.USub,ast.UAdd,ast.Constant,ast.Name,ast.Load)
    assert all(isinstance(t,permitted) for t in ast.walk(tree))
    return compile(ast.fix_missing_locations(ExactDivision().visit(tree)),'<coefficient>','eval')

def projection(poly):
    """Exact monomial moments in the plain triple Laguerre basis."""
    out={}
    for p,c in S.Poly(poly,*O.xyz).terms():
        c=Q(int(S.numer(c)),int(S.denom(c)))
        import itertools
        for j in itertools.product(*(range(k+1) for k in p)):
            factor=math.prod((-1)**idx*math.factorial(power)*math.comb(power,idx)
                             for power,idx in zip(p,j))
            out[j]=out.get(j,Q(0))+c*factor
    return {j:c for j,c in out.items() if c}

def main():
    minimal='--minimal' in sys.argv
    raw=json.loads((OUT/('minimal_sigma_corrections.json' if minimal else 'codex_recurrence.json')).read_text())
    differential=json.loads((OUT/('minimal_sigma_operators.json' if minimal else 'cleared_operators.json')).read_text())
    rng=random.Random(20261021)
    sources=[((0,0,0),Q(rng.randint(1,5),3)),
             ((1,2,0),Q(rng.randint(1,5),7)),
             ((2,0,1),Q(rng.randint(1,5),11))]
    poly=sum(S.Rational(c.numerator,c.denominator)*S.prod(S.laguerre(n,z) for n,z in zip(j,O.xyz))
             for j,c in sources)
    passed=0
    # Pi is a formal multiplicative constant in the affected blocks. Substitution
    # by 22/7 makes these exact rational tests; direct polylog tests use true pi.
    for kind,blocks in raw.items():
        for key,stencil in blocks.items():
            op={ast.literal_eval(d):R.parse(f) for d,f in differential[kind][key].items()}
            par={O.eps:S.Rational(3,2),O.E:S.Rational(-7,3),O.Z:2,O.lam:S.Rational(4,5),S.pi:S.Rational(22,7)}
            expected=projection(sum(f.subs(par)*S.diff(poly,*sum(([z,n] for z,n in zip(O.xyz,d)),[]))
                                    for d,f in op.items()))
            actual={}
            constants={'eps':Q(3,2),'E':Q(-7,3),'Z':Q(2),'lam':Q(4,5),'pi':Q(22,7),'Q':Q}
            for offset,expression in stencil.items():
                offset=ast.literal_eval(offset)
                targets=[(tuple(j[i]-offset[i] for i in range(3)),c) for j,c in sources]
                targets=[(j,c) for j,c in targets if min(j)>=0]
                if not targets:continue
                code=compile_coefficient(expression)
                for target,c in targets:
                    env=dict(constants,l=target[0],m=target[1],n=target[2])
                    value=eval(code,{'__builtins__':{}},env)*c
                    actual[target]=actual.get(target,Q(0))+value
            actual={j:c for j,c in actual.items() if c}
            assert actual==expected,(kind,key,[(j,actual.get(j,0)-expected.get(j,0))
                       for j in set(actual)|set(expected) if actual.get(j,0)!=expected.get(j,0)][:5])
            passed+=1
            print('PASS',kind,key,len(expected),'output Laguerre coefficients',flush=True)
    print('PASS: all',passed,'saved block stencils, exact random rational finite combination')

if __name__=='__main__':main()
