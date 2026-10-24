"""Create the report's fit tables directly from the saved serial outputs."""
from pathlib import Path
O=Path(__file__).resolve().parent
fp=[l.split() for l in (O/'x19_fit40_results.txt').read_text().splitlines() if l and not l.startswith('#')]
bp=[l.split() for l in (O/'x19_fit40_big_results.txt').read_text().splitlines() if l and not l.startswith('#')]
f=lambda x:f'{float(x):.3e}'
text='''**VERIFIED numerical tests, not a nonmembership theorem.** `spec40.txt`
contains C0,...,C4 at E=a21=0. The fitted targets are exactly its C1,C2,C3.
All five order-two atoms, all eighteen classical generators and all 4*23 x16
states were evaluated independently in `x19_atoms.jl` from their defining
functions/integrals. The full family has 115 columns before polynomial
multipliers; duplicates are retained and handled by rank-revealing factorizations.
No symmetry restriction is imposed on the fitting span.

On R=1 use all monomials r1^i r2^j r12^k with j=0 or 1 and i+j+k<=d.
This is the complete degree-at-most-d polynomial restriction, since
r2²=1-r1², and contains (d+1)² monomials. These are slice membership tests;
by themselves they do not assert a homogeneous global distance representation.

The fixed geometric cut is 0.25<alpha,theta<pi-0.25, selected without looking
at target values. Of the 4,000 supplied points it retains 2,801: 2,240 training
and 561 held-out. Float64 panel rules with 256 and 512 nodes agree over every
retained point to max scaled discrepancy 1.034e-14. The earlier 128-node rule
failed the guard; its values were never fitted. Its worst state's convergence
sequence is saved in `x19_quadrature.out`.

The following table gives Float64 **held-out relative L2 residuals**, with an
SVD cutoff 1e-12 relative to the largest singular value. Numerical rank is a
conditioning diagnostic, not a rational-function rank certificate.

| Family | Degree | Columns | Numerical rank | C1 | C2 | C3 |
|---|---:|---:|---:|---:|---:|---:|
'''
for family in ('order2','order2_classical','all115'):
 for d in sorted(set(int(r[1]) for r in fp if r[0]==family)):
  rr=sorted([r for r in fp if r[0]==family and int(r[1])==d],key=lambda r:int(r[5]))
  text+='| '+' | '.join([family,str(d),rr[0][2],rr[0][3]]+[f(r[7]) for r in rr])+' |\n'
text+='''
The SVD cutoff sensitivity at d=2 is material: 1e-10 / 1e-12 / 1e-14
produce different effective ranks and residuals. Even the known degree-four
psi41 control has held-out residual 1.833e-11 in Float64 at the main cutoff.
Thus residuals near 1e-11 in that table cannot distinguish an identity.
Every training residual, held-out relative maximum, cutoff and rank is retained
in `x19_fit40_results.txt`.

**256-bit follow-up.** A deterministic 1,000-point subset has 800 training and
200 separate held-out points. Each panel uses a Gauss rule on panels graded to
the arctan poles and the complex quartic-root distance. 96 vs 128 nodes per
panel agree to 1.906e-58 at five validation points. A second implementation,
adaptive mpmath quadrature of the prior independent x16 actual-integral
checker, compares all 92 states at two fresh cache points; see
`x19_panelcheck.json` for its errors. This guards against a shared finite-rule
convergence error.

Column-scaled pivoted BigFloat QR uses a relative R-diagonal threshold 1e-42.
The table shows **training / held-out relative L2 residuals**. The full degree-1
family is still ill-conditioned and extrapolates substantially; its small
training residual must not be treated as an identity.

| Family | Degree | Columns | QR numerical rank | C1 train / held | C2 train / held | C3 train / held |
|---|---:|---:|---:|---:|---:|---:|
'''
for family in ('all115_big','order2_big'):
 for d in sorted(set(int(r[1]) for r in bp if r[0]==family)):
  rr=sorted([r for r in bp if r[0]==family and int(r[1])==d and int(r[4])<=3],key=lambda r:int(r[4]))
  text+='| '+' | '.join([family,str(d),rr[0][2],rr[0][3]]+[f(r[5])+' / '+f(r[6]) for r in rr])+' |\n'
text+='''
The degree-four psi41 control has full QR rank 125 and residual
1.720e-75 training / 1.974e-75 held-out. This confirms that the high-precision
pipeline can recover a known identity at the working precision; the psi40
residuals in the tested spans are much larger. The initial unrestricted
high-precision sweep stopped at degree one (460 columns). A symmetry reduction
then enabled the full singlet degree-two search below. Higher degrees in the
full family have only the Float64 evidence, where rank truncation is severe.

SINGLET_TABLE_PLACEHOLDER

**Denominator controls.** Put

    delta=4 r1² r2²-(1-r12²)²,
    Dclear=r1³ r2³ r12³ (2-r12²)² delta²  (R=1).

The homogeneous version includes an extra factor R². The exact geometry checker
verifies that it is divisible by every one of the 92 exported first-order x16
row clearings. We fit delta*Cj, delta²*Cj and Dclear*Cj as well as Cj.
The following are held-out relative L2 residuals in the same 115-state span:

| Multiplier | Precision / degree | C1 | C2 | C3 |
|---|---|---:|---:|---:|
'''
for label in ('all115_delta','all115_delta2','all115_commonclear'):
 rr=sorted([r for r in fp if r[0]==label and r[1]=='3'],key=lambda r:int(r[5]))
 text+='| '+' | '.join([label,'Float64 / 3']+[f(r[7]) for r in rr])+' |\n'
for rhs0,label in [(4,'delta'),(7,'delta2'),(10,'Dclear')]:
 rr=sorted([r for r in bp if r[0]=='all115_big' and r[1]=='1' and rhs0<=int(r[4])<rhs0+3],key=lambda r:int(r[4]))
 text+='| '+' | '.join([label,'256 bit / 1']+[f(r[6]) for r in rr])+' |\n'
text+='''
These are residuals for the multiplied targets, hence weighted residuals for
the original functions. Multiplication raises the required numerator degree;
failure at these low degrees does not exclude a rational representation with
these denominators. The existing derivative clearings need not be the only
possible denominators of a membership relation. The tests specifically try the
requested discriminant/row-clear options and **do not** conflate their failure
with a proof that every denominator choice fails.

**Assessment:** no exact membership of C1,C2,C3 in the tested old spans was found.
Their approximation improves quickly with degree. Neither an exact old-module
representation nor its impossibility has been established. The high-precision
control shows that these residuals are not merely the arithmetic floor. Exact
nonmembership even at a tested finite degree is not certified: discarded tiny
directions could matter with large coefficients. Higher degrees and other
denominator choices remain possible. No coefficient
PSLQ reconstruction or claimed closed form is based on these approximate fits.

'''
sp=[l.split() for l in (O/'x19_fit40_singlet_results.txt').read_text().splitlines() if l and not l.startswith('#')]
st="""**Full degree-two singlet follow-up at 256 bits.** The atom family is closed
under electron exchange with constant coefficients. Center alpha by pi/2; the
order-two atoms have exchange parities (+,+,-,+,-). Classical L/T pairs
exchange with a minus sign, the two remaining L/log pairs with a plus sign;
classical a is odd and p,a²,p²,1 are even. The panel families exchange
1<->2 and 3<->4. Combining paired sums/differences gives 60 even and 55 odd
generators, including presentation dependencies. Direct actual-function checks
at two points verify these parities to 9.0e-77.

Let s=r1+r2 and d=r1-r2. Since s²+d²=2 on R=1, every degree-at-most-two
polynomial reduces to s^i d^j r12^k with j=0 or 1. The complete singlet span is
thus even atoms times six even monomials plus odd atoms times three odd
monomials: 60*6+55*3=525 columns. This is the singlet image of the full
1,035-column degree-two family, not a chosen subset of its atoms. The same
800/200 point split and QR cutoff give numerical rank 463.

| Multiplier | C1 train / held | C2 train / held | C3 train / held |
|---|---:|---:|---:|
"""
for lo,lab in [(1,'1'),(4,'delta'),(7,'delta2'),(10,'Dclear')]:
 rr=[r for r in sp if lo<=int(r[3])<lo+3]
 st+='| '+' | '.join([lab]+[f(r[4])+' / '+f(r[5]) for r in rr])+' |\n'
st+="""
A separately manufactured member of this full panel-polynomial span has
training residual 1.328e-47 and held-out residual 2.425e-38. The fit remains
ill-conditioned and amplifies errors on held-out points, but the Cj residuals
are far above this control. The complete output is in
`x19_fit40_singlet_results.txt`; peak RSS was 806,502,400 bytes. Full-family
high-precision singlet degree>=3, unrestricted exact rational-rank certification
and unbounded-degree membership remain **OPEN**.
"""
text=text.replace('SINGLET_TABLE_PLACEHOLDER',st)
p=O/'x19_report.md';s=p.read_text();start=s.index('## 2. ψ40 charge-sector membership tests');end=s.index('## 3. Line charges',start)
s=s[:start]+'## 2. ψ40 charge-sector membership tests\n\n'+text+'\n'+s[end:];p.write_text(s)
