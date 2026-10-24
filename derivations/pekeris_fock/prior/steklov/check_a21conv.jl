include("steklov.jl")
model = channel_model(6)
mom = gauge_moments(model, 2)
println("2<U^2/2,Y2> quadrature = ", mom[3][2, 1], "   formula = ", 2 * (pi + 4) / (6pi) - 1 / 8)
for nq in (40, 80, 120)
    m2 = channel_model(6; nquad = nq)
    println("nquad=$nq: |dU|=", norm(m2.U - model.U, Inf), " |dW|=", norm(m2.W - model.W, Inf), " |dC|=", norm(m2.C - model.C, Inf))
end
