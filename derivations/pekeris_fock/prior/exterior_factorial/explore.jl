include("/Users/nathompson7/Documents/Codex/2026-09-05/how-much-compu/outputs/boundary_matching/cusp_matching.jl")
for d in (2, 4, 6)
    b = cuspmodel(d)
    z = cuspsolve(b)
    m = length(b.K)
    T = b.W + 2z.E * I
    ev = eigvals([zeros(m, m) I; -T -2b.U])
    println("d=$d m=$m E=$(z.E) symU=$(norm(b.U-b.U')) symT=$(norm(T-T')) symD=$(norm(b.D-b.D'))")
    println(sort(ev, by = real))
end
