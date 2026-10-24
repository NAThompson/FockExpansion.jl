include("atoms_ritz.jl")
setprecision(BigFloat, 320)
nd = nodes(BigFloat(1) / 16, 72)
v, gx, gy = atom_table(psi30h_xy, nd; float64 = true)
bad = [q for q in eachindex(v) if !isfinite(v[q]) || !isfinite(gx[q]) || !isfinite(gy[q])]
println(length(bad), " bad of ", length(v))
for q in bad[1:min(8, end)]
    α, μ = angles(nd.x[q], nd.y[q])
    println("x=", Float64(nd.x[q]), " y=", Float64(nd.y[q]), " α=", Float64(α), " θ=", Float64(theta(nd.x[q], nd.y[q])), " w=", Float64(nd.w[q]), " v=", Float64(v[q]), " gx=", Float64(gx[q]))
end
println("max weight among bad: ", isempty(bad) ? 0 : Float64(maximum(nd.w[bad])))
