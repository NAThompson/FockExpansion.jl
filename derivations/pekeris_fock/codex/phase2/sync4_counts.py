import sys
sys.dont_write_bytecode=True
import json,ast
from pathlib import Path
P=Path(__file__).resolve().parent
c=json.loads((P.parents[1]/'claude2/order2_recurrence.json').read_text())
o=json.loads((P/'counts.json').read_text())
mapname={'Om':'U','al':'x','la':'L','lR':'Slog','one':'one','beP':'V','P':'V','lb':'U'}
lines=['# Saved eight-atom counts, independently parsed','', '| Operator | Target <- source | Offsets | Shell range | Related six-atom block/count |','|---|---|---:|---|---|']
for kind in ('sigma2P','C1','C2'):
 for key,st in c[kind].items():
  target,source=ast.literal_eval(key);sh=[sum(ast.literal_eval(k)) for k in st]
  key2='('+mapname[target]+','+mapname[source]+')'
  other=o['S2P' if kind=='sigma2P' else kind].get(key2)
  # These are related blocks, not generally equal operators.
  related=key2+' / '+str(other['fixed_count']) if other else '—'
  lines.append(f'| {kind} | {target} <- {source} | {len(st)} | [{min(sh)}, {max(sh)}] | {related} |')
lines+=['','Related means only basis correspondence; U=Om+pi lb and V=(beP+pi P/2)/r12 mix blocks and derivatives. Compare after removing row factors.','', 'Row factors: `'+str(c['D'])+'`.']
(P/'sync4_counts.md').write_text('\n'.join(lines)+'\n')
print({k:len(c[k]) for k in ('sigma2P','C1','C2')})
print('total',sum(len(st) for k in ('sigma2P','C1','C2') for st in c[k].values()))
print('old stats example',next(iter(o['S2P'].values())))
