include("steklov.jl")
module CodexRef
include("/Users/nathompson7/Documents/Codex/2026-09-05/how-much-compu/outputs/boundary_matching/cusp_matching.jl")
end

for d in (2, 6)
    rule = sphere_rule(40, 40)
    basis = harmonic_basis(d)
    A = harmonic_values(basis, sphere_q(rule), sphere_s(rule))
    println("d=$d sum w = ", sum(rule.w), "  gram err = ", norm(A' * (rule.w .* A) - I, Inf))
    mine = channel_model(d)
    ref = CodexRef.cuspmodel(d; nquad = 100)
    println("  U diff ", norm(mine.U - ref.U, Inf), "  D diff ", norm(mine.D - ref.D, Inf),
            "  W diff ", norm(mine.W - ref.W, Inf), "  derr ", mine.derr)
    Vref = CodexRef.angular(d).V
    println("  C diff ", norm(mine.C - Vref, Inf))
end
