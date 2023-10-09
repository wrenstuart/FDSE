using CairoMakie

#=
f = Figure(backgroundcolor = :white, resolution=(800,300))
ax = Axis(f[1,1])
=#

x = range(0,10,length=100)
y=sin.(x)
scatter(x,y)
#f, ax, scat = scatter(x,y)
#f
