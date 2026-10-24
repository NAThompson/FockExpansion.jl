include("fitbig.jl")
pts, vals = load("spec41.txt")
pts, vals = pts[1:1600], vals[1:1600]
names = ["1", "α", "β", "l1", "l2", "l3", "Ω", "P", "Pα", "Pβ", "Pl1", "Pl2", "Pl3", "PΩ"]
full = collect(1:14)
for drop in full
    sel = setdiff(full, [drop])
    print("drop ", names[drop], ": ")
    runfit(pts, vals, sel, 6)
end
