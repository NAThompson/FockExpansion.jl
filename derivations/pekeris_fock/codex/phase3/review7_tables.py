import pathlib,re,math,json,hashlib
root=pathlib.Path(__file__).resolve().parents[2]
paths=[root/'claude3/atoms3.jl',root/'claude3/run3.jl',root/'claude3/check3.out',root/'claude3/param_check.out',*sorted((root/'claude3').glob('run3_*.out')),*sorted((root/'claude3').glob('a21_*.out')),root/'claude3/states/min5_N4.txt',root/'claude3/states/min5_N8.txt']
manifest={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
(root/'codex/phase3/review7_input_hashes.json').write_text(json.dumps(manifest,indent=2)+'\n')
rows=[]
for line in (root/'claude3/run3_min5.out').read_text().splitlines():
 m=re.search(r'N=\s*(\d+) size=\s*(\d+)\s+E-E0=\s*(\S+)',line)
 if m: rows.append((int(m[1]),int(m[2]),float(m[3])))
print('N interval | log(error)/log(size) slope | log(error)/size slope')
for a,b in zip(rows,rows[1:]):
 d=math.log(b[2]/a[2]);print(a[0],b[0],d/math.log(b[1]/a[1]),d/(b[1]-a[1]))
print('Semilog convex iff successive last-column secants increase.')
ref=.47674787900
for name in ('a21_min5.out','a21_a2psi.out','a21_min5psi.out'):
 print(name)
 for line in (root/'claude3'/name).read_text().splitlines():
  a=float(re.search(r'a21_global=(\S+)',line)[1]);print(line.split()[0], 'reference difference',a-ref)
