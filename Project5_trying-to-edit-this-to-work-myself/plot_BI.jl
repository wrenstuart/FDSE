
# This script reads in output from KH.jl, makes a plot, and saves an animation

using Oceananigans, Oceananigans.Units, JLD2, Plots, Printf
Lx = 100kilometers
Ly = 300kilometers

# Set the filename (without the extension)
filename = "BI"

# Read in the first iteration.  We do this to load the grid
# filename * ".jld2" concatenates the extension to the end of the filename
ζ_ic = FieldTimeSeries(filename * ".jld2", "ζ", iterations = 0)
S_ic = FieldTimeSeries(filename * ".jld2", "S", iterations = 0)
#=u_ic = FieldTimeSeries(filename * ".jld2", "u", iterations = 0)
v_ic = FieldTimeSeries(filename * ".jld2", "v", iterations = 0)
w_ic = FieldTimeSeries(filename * ".jld2", "w", iterations = 0)
b_ic = FieldTimeSeries(filename * ".jld2", "b", iterations = 0)
ω_ic = FieldTimeSeries(filename * ".jld2", "ω", iterations = 0)
χ_ic = FieldTimeSeries(filename * ".jld2", "χ", iterations = 0)
ϵ_ic = FieldTimeSeries(filename * ".jld2", "ϵ", iterations = 0)=#

## Load in coordinate arrays
## We do this separately for each variable since Oceananigans uses a staggered grid
xζ, yζ, zζ = nodes(ζ_ic)
xS, yS, zS = nodes(S_ic)
#=xu, yu, zu = nodes(u_ic)
xv, yv, zv = nodes(v_ic)
xw, yw, zw = nodes(w_ic)
xb, yb, zb = nodes(b_ic)
xω, yω, zω = nodes(ω_ic)
xχ, yχ, zχ = nodes(χ_ic)
xϵ, yϵ, zϵ = nodes(ϵ_ic)=#

## Now, open the file with our data
file_xy = jldopen(filename * ".jld2")

## Extract a vector of iterations
iterations = parse.(Int, keys(file_xy["timeseries/t"]))

@info "Making an animation from saved data..."

# Here, we loop over all iterations
anim = @animate for (i, iter) in enumerate(iterations)

    @info "Drawing frame $i from iteration $iter..."

    ζ_xz = file_xy["timeseries/ζ/$iter"][:, :, 1];
    S_xz = file_xy["timeseries/S/$iter"][:, :, 1];
    #=u_xz = file_xz["timeseries/u/$iter"][:, 1, :];
    v_xz = file_xz["timeseries/v/$iter"][:, 1, :];
    w_xz = file_xz["timeseries/w/$iter"][:, 1, :];
    b_xz = file_xz["timeseries/b/$iter"][:, 1, :];
    ω_xz = file_xz["timeseries/ω/$iter"][:, 1, :];
    χ_xz = file_xz["timeseries/χ/$iter"][:, 1, :];
    ϵ_xz = file_xz["timeseries/ϵ/$iter"][:, 1, :];=#

    t = file_xy["timeseries/t/$iter"];

        ζ_xy_plot = heatmap(xζ, yζ, ζ_xz'; color = :thermal, xlabel = "x", ylabel = "y", aspect_ratio = :equal, xlims = (0, Lx), ylims = (0, Ly),clims=(-3e-4,3e-4)); 
        S_xy_plot = heatmap(xS, yS, S_xz'; color = :thermal, xlabel = "x", ylabel = "y", aspect_ratio = :equal, xlims = (0, Lx), ylims = (0, Ly)); 
        #=b_xz_plot = heatmap(xb, zb, b_xz'; color = :thermal, xlabel = "x", ylabel = "z", aspect_ratio = :equal, xlims = (0, Lx), ylims = (0, Lz)); 
        ω_xz_plot = heatmap(xω, zω, ω_xz'; color = :thermal, xlabel = "x", ylabel = "z", aspect_ratio = :equal, xlims = (0, Lx), ylims = (0, Lz)); 
        χ_xz_plot = heatmap(xχ, zχ, χ_xz'; color = :thermal, xlabel = "x", ylabel = "z", aspect_ratio = :equal, xlims = (0, Lx), ylims = (0, Lz)); 
        ϵ_xz_plot = heatmap(xϵ, zϵ, ϵ_xz'; color = :thermal, xlabel = "x", ylabel = "z", aspect_ratio = :equal, xlims = (0, Lx), ylims = (0, Lz)); =#

    ζ_title = @sprintf("ζ, t = %s", round(t));
    S_title = @sprintf("S, t = %s", round(t));
    #=u_title = @sprintf("u, t = %s", round(t));
    v_title = @sprintf("v, t = %s", round(t));
    w_title = @sprintf("w, t = %s", round(t));
    b_title = @sprintf("b, t = %s", round(t));
    ω_title = @sprintf("vorticity (ω), t = %s", round(t));
    ϵ_title = @sprintf("KE dissipation (ϵ), t = %s", round(t));
    χ_title = @sprintf("buoyancy variance dissipation (χ), t = %s", round(t));=#

# Combine the sub-plots into a single figure
    #plot(b_xz_plot, ω_xz_plot, ϵ_xz_plot, χ_xz_plot, layout = (4, 1), size = (1200, 800),
    plot(ζ_xy_plot, S_xy_plot, layout = (1,2))
    #title = [b_title ω_title ϵ_title χ_title])
    title = [ζ_title S_title]

    iter == iterations[end] && close(file_xy)
end

close(file_xy)

# Save the animation to a file
mp4(anim, "BI.mp4", fps = 20) # hide
