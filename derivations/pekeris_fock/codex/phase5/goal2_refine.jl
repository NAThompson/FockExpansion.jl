include("goal2_panel_rank.jl")
for p=((.3,-.6),(-.27,-.7)); experiment(p...;D=10,M=40,rho=.12); end
