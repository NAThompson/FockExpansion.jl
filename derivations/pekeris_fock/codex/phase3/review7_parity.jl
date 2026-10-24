include("review7_independent.jl")
for (x,y) in ((big"0.2",big"0.6"),(big"0.65",big"0.8"))
 for g in (:Lt,:alp,:V,:U)
  sign=g in (:alp,:U) ? -1 : 1
  err=abs(atom(g,-x,y)-sign*atom(g,x,y))
  @printf("x=%.2f y=%.2f atom=%s exchange residual=%.3e\n",Float64(x),Float64(y),g,Float64(err))
  @assert err<big"1e-65"
 end
end
