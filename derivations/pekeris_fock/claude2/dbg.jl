include("atoms_ritz.jl")
setprecision(BigFloat, 256)
nd = nodes(BigFloat(1) / 12, 50)
for q in eachindex(nd.x)
    try
        ForwardDiff.gradient(z -> chi_xy(z[1], z[2]), [nd.x[q], nd.y[q]])
    catch e
        println("fail at x=", Float64(nd.x[q]), " y=", Float64(nd.y[q]), " y-x=", Float64(nd.y[q] - nd.x[q]), " : ", sprint(showerror, e)[1:min(200, end)])
        break
    end
end
