using CairoMakie

#=
f = Figure(backgroundcolor = :white, resolution=(800,300))
ax = Axis(f[1,1])
=#

x = range(0,10,length=100)
y = range(0,10,length=100)
out_x = cos.(x)
out_y = sin.(y)
out_xy = out_x*out_y'

#heatmap(x,y,out_xy)

f = Figure(backgroundcolor=:white, reolsution=(800,300))
ax = Axis(f[1,1])
heatmap!(ax,x,y,out_xy)
f