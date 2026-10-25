ps = [2147483647, 2147483629, 2147483587]
files = {2147483647: 'jets_p1.txt', 2147483629: 'jets_2147483629.txt', 2147483587: 'jets_2147483587.txt'}
data = {p: {tuple(map(int, l.split()[:2])): [ZZ(x) for x in l.split()[2:]] for l in open(files[p])} for p in ps}
out = open('jets_Q.txt', 'w')
for key in sorted(data[ps[0]]):
    vals = [rational_reconstruction(crt([data[p][key][i] for p in ps], ps), prod(ps)) for i in range(5)]
    out.write('%d %d %s\n' % (key[0], key[1], ' '.join(str(v) for v in vals)))
    print(key, [len(str(v)) for v in vals])
