include("/Users/nathompson7/Documents/Codex/2026-09-05/how-much-compu/outputs/boundary_matching/cusp_matching.jl")
using LinearAlgebra
E = -2.903724377034
for d in 2:2:16
    b = cuspmodel(d)
    m = length(b.K)
    T = b.W + 2E * I
    ev = sort(real(eigvals([zeros(m, m) I; -T -2b.U])))
    println("d=$d m=$m slow=[", round(ev[1]; digits=3), ", ", round(ev[m]; digits=3),
            "] fast_min=", round(ev[m+1]; digits=3), " gap=", round(ev[m+1] - ev[m]; digits=3))
end
