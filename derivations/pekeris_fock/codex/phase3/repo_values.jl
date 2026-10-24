using FockExpansion
setprecision(BigFloat,192)
for (aa,tt,zz,ee,cc) in (("0.9","1.1","2","-2.9","0.123"),("2.0","0.7","1.7","-2.2","0.31"),("1.2","1.8","0","-1.3","0.27"))
 a,t,z,e,c=parse.(BigFloat,(aa,tt,zz,ee,cc))
 p=psi30_parts(a,t;Z=z,E=e,a21=c,rtol=big"1e-45")
 pf=psi30_parts(Float64(a),Float64(t);Z=Float64(z),E=Float64(e),a21=Float64(c),rtol=1e-12)
 println(join((aa,tt,zz,ee,cc,string(p.classical),string(p.state),string(p.elliptic),string(p.value),string(pf.value)),"\t"));flush(stdout)
end
