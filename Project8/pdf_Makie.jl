# This script reads in output from KH.jl, makes a plot, and saves an animation

using Oceananigans, JLD2, Printf, CairoMakie

CairoMakie.activate!()

pre = "Project8/raw-output/BI_xz"
Ris = [0, 1, 2, 3, 4]
post = ""

for i = 1 : 5

log_Ri = Ris[i]
Ri = 10 ^ log_Ri
Ro = (1+Ri) ^ -0.5
lim = 10 * Ro

α = 2

# Set the filename (without the extension)
filename = pre * "_Ri_" * string(log_Ri) * "_α_" * string(α)

# Read in the first iteration.  We do this to load the grid
# filename * ".jld2" concatenates the extension to the end of the filename
u_ic = FieldTimeSeries(filename * ".jld2", "u", iterations = 0)
v_ic = FieldTimeSeries(filename * ".jld2", "v", iterations = 0)
w_ic = FieldTimeSeries(filename * ".jld2", "w", iterations = 0)
ζ_ic = FieldTimeSeries(filename * ".jld2", "ζ", iterations = 0)
δ_ic = FieldTimeSeries(filename * ".jld2", "δ", iterations = 0)

## Load in coordinate arrays
## We do this separately for each variable since Oceananigans uses a staggered grid
xu, yu, zu = nodes(u_ic)
xv, yv, zv = nodes(v_ic)
xw, yw, zw = nodes(w_ic)
xζ, yζ, zζ = nodes(ζ_ic)
xδ, yδ, zδ = nodes(δ_ic)

## Now, open the file with our data
file_xy = jldopen(filename * ".jld2")

## Extract a vector of iterations
iterations = parse.(Int, keys(file_xy["timeseries/t"]))

@info "Making an animation from saved data..."
@info length(iterations)/20
@info iterations[end]

iter = iterations[length(iterations)]
ζ_xy = file_xy["timeseries/ζ/$iter"][:, :, 1]

f = Figure(backgroundcolor=:white, reolsution=(800,300))
ax = Axis(f[1,1], limits = ((-3, 3), (0, 5)))
#=density!(ax,vec(ζ_xy),bins=100)
ax = Axis(f[1,2])
heatmap!(ax,xζ,yζ,ζ_xy,colormap=:balance)
=#


#=
function update_frame(iter)
    i, time = iter;
    @info i
    ζ_xy = file_xy["timeseries/ζ/$time"][:, :, 1]*1e4;
    #ζ_xy = cos.(xζ.+Float64(i)/20)*(sin.(yζ))'
    @info findmax(ζ_xy)[1]
    hm = heatmap!(ax,xζ,yζ,ζ_xy,colormap=:balance,colorrange=(-3,3))
    Colorbar(f[1,2],hm)
end
=#

function update_frame(iter)
    i, time = iter;
    @info i
    empty!(ax)
    ζ_xy = file_xy["timeseries/ζ/$time"][:, :, 1] * 1e4 * Ri;
    CairoMakie.density!(ax, vec(ζ_xy), color = :tomato)
end

record(update_frame, f, "Project8/" * "hist_Ri_" * string(log_Ri) * "_α_" * string(α) * ".mp4", enumerate(iterations[1:1:end]); framerate=20)

close(file_xy)
nothing

end