# This script reads in output from KH.jl, makes a plot, and saves an animation

using Oceananigans, JLD2, Printf, CairoMakie

CairoMakie.activate!()

pre = "Project7/raw-output/BI_xz"
Ris = [0, 1, 2, 3, 4]

for i = 1 : 5

log_Ri = Ris[i]
Ri = 10 ^ log_Ri
Ro = (1+Ri) ^ -0.5
lim = 10 * Ro

if i == 5
    λ = 0
else
    λ = 2
end

# Set the filename (without the extension)
filename = pre * "_Ri-" * string(log_Ri) * "_λ-" * string(λ)

# Read in the first iteration.  We do this to load the grid
# filename * ".jld2" concatenates the extension to the end of the filename
u_ic = FieldTimeSeries(filename * ".jld2", "u", iterations = 0)
v_ic = FieldTimeSeries(filename * ".jld2", "v", iterations = 0)
w_ic = FieldTimeSeries(filename * ".jld2", "w", iterations = 0)
#T_ic = FieldTimeSeries(filename * ".jld2", "T", iterations = 0)
#S_ic = FieldTimeSeries(filename * ".jld2", "S", iterations = 0)
ζ_ic = FieldTimeSeries(filename * ".jld2", "ζ", iterations = 0)

## Load in coordinate arrays
## We do this separately for each variable since Oceananigans uses a staggered grid
xu, yu, zu = nodes(u_ic)
xv, yv, zv = nodes(v_ic)
xw, yw, zw = nodes(w_ic)
#xT, yT, zT = nodes(T_ic)
#xS, yS, zS = nodes(S_ic)
xζ, yζ, zζ = nodes(ζ_ic)

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
ax = Axis(f[1,1], limits = ((-10, 10), (0, 1)))
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
    ζ_xy = file_xy["timeseries/ζ/$time"][:, :, 1] * 1e4 / Ro;
    CairoMakie.density!(ax, vec(ζ_xy), color = :tomato)
end

record(update_frame, f, "Project7/" * "hist_Ri-" * string(log_Ri) * "_λ-" * string(λ) * ".mp4", enumerate(iterations[1:1:end]); framerate=20)

#=
# Here, we loop over all iterations
anim = @animate for (i, iter) in enumerate(iterations)

    @info "Drawing frame $i from iteration $iter..."

    u_xz = file_xz["timeseries/u/$iter"][:, 1, :];
    v_xz = file_xz["timeseries/v/$iter"][:, 1, :];
    w_xz = file_xz["timeseries/w/$iter"][:, 1, :];
    b_xz = file_xz["timeseries/b/$iter"][:, 1, :];
    ω_xz = file_xz["timeseries/ω/$iter"][:, 1, :];
    χ_xz = file_xz["timeseries/χ/$iter"][:, 1, :];
    ϵ_xz = file_xz["timeseries/ϵ/$iter"][:, 1, :];

    t = file_xz["timeseries/t/$iter"];

        b_xz_plot = heatmap(xb, zb, b_xz'; color = :thermal, xlabel = "x", ylabel = "z", aspect_ratio = :equal, xlims = (0, Lx), ylims = (0, Lz)); 
        ω_xz_plot = heatmap(xω, zω, ω_xz'; color = :thermal, xlabel = "x", ylabel = "z", aspect_ratio = :equal, xlims = (0, Lx), ylims = (0, Lz)); 
        χ_xz_plot = heatmap(xχ, zχ, χ_xz'; color = :thermal, xlabel = "x", ylabel = "z", aspect_ratio = :equal, xlims = (0, Lx), ylims = (0, Lz)); 
        ϵ_xz_plot = heatmap(xϵ, zϵ, ϵ_xz'; color = :thermal, xlabel = "x", ylabel = "z", aspect_ratio = :equal, xlims = (0, Lx), ylims = (0, Lz)); 

    u_title = @sprintf("u, t = %s", round(t));
    v_title = @sprintf("v, t = %s", round(t));
    w_title = @sprintf("w, t = %s", round(t));
    b_title = @sprintf("b, t = %s", round(t));
    ω_title = @sprintf("vorticity (ω), t = %s", round(t));
    ϵ_title = @sprintf("KE dissipation (ϵ), t = %s", round(t));
    χ_title = @sprintf("buoyancy variance dissipation (χ), t = %s", round(t));

# Combine the sub-plots into a single figure
    plot(b_xz_plot, ω_xz_plot, ϵ_xz_plot, χ_xz_plot, layout = (4, 1), size = (1200, 800),
    title = [b_title ω_title ϵ_title χ_title])

    iter == iterations[end] && close(file_xz)
end

=#
close(file_xy)
nothing

end