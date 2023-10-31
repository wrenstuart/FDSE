
# This script reads in output from gravitycurrent.jl, makes a plot, and saves an animation

using Oceananigans, Oceananigans.Units, JLD2, Plots, Printf
# given f, simnum


# Set the filename (without the extension)
filename_xz = "Project6/BI_xz_$simnum"
filename_xy = "Project6/BI_xy_$simnum"

# Read in the first iteration.  We do this to load the grid
# filename * ".jld2" concatenates the extension to the end of the filename
u_ic = FieldTimeSeries(filename_xy * ".jld2", "u", iterations = 0)
v_ic = FieldTimeSeries(filename_xy * ".jld2", "v", iterations = 0)
w_ic = FieldTimeSeries(filename_xy * ".jld2", "w", iterations = 0)
T_ic = FieldTimeSeries(filename_xy * ".jld2", "T", iterations = 0)
S_ic = FieldTimeSeries(filename_xy * ".jld2", "S", iterations = 0)
ζ_ic = FieldTimeSeries(filename_xy * ".jld2", "ζ", iterations = 0)

## Load in coordinate arrays
## We do this separately for each variable since Oceananigans uses a staggered grid
xu, yu, zu = nodes(u_ic)
xv, yv, zv = nodes(v_ic)
xw, yw, zw = nodes(w_ic)
xT, yT, zT = nodes(T_ic)
xS, yS, zS = nodes(S_ic)
xζ, yζ, zζ = nodes(ζ_ic)

## Now, open the file with our data
file_xy = jldopen(filename_xy * ".jld2")
file_xz = jldopen(filename_xz * ".jld2")

## Extract a vector of iterations
iterations = parse.(Int, keys(file_xy["timeseries/t"]))

@info "Making an animation from saved data..."

t_save = zeros(length(iterations))

# Here, we loop over all iterations
anim = @animate for (i, iter) in enumerate(iterations)

    @info "Drawing frame $i from iteration $iter..."

    #T_xz = file_xz["timeseries/T/$iter"][:, 1, :];
    #S_xz = file_xz["timeseries/S/$iter"][:, 1, :];

# If you want an x-y slice, you can get it this way:
    T_xy = file_xy["timeseries/T/$iter"][:, :, 1];
    S_xy = file_xy["timeseries/S/$iter"][:, :, 1];
    ζ_xy = file_xy["timeseries/ζ/$iter"][:, :, 1];

    t = file_xy["timeseries/t/$iter"];

    # Save some variables to plot at the end
    t_save[i] = t # save the time

        T_xy_plot = heatmap(xT/1kilometer, yT/1kilometer, T_xy'; color = :thermal, xlabel = "x (km)", ylabel = "y (km)")#, aspect_ratio = :equal);  
        S_xy_plot = heatmap(xS/1kilometer, yS/1kilometer, S_xy'; color = :haline, xlabel = "x (km)", ylabel = "y (km)")#, aspect_ratio = :equal);  
        ζ_xy_plot = heatmap(xζ/1kilometer, yζ/1kilometer, ζ_xy'/f; color = :balance, xlabel = "x (km)", ylabel = "y (km)")#, aspect_ratio = :equal, clims=(-1.5,1.5));  

    T_title = @sprintf("Temperature");
    S_title = @sprintf("Salinity");
    ζ_title = @sprintf("ζ/f");

# Combine the sub-plots into a single figure
    plot(T_xy_plot, S_xy_plot, ζ_xy_plot, layout = (1, 3), size = (802, 600),
    title = [T_title S_title ζ_title])

    iter == iterations[end] && close(file_xy)
end

close(file_xz)
close(file_xy)

# Save the animation to a file
mp4(anim, "Project6/BI_$simnum.mp4", fps = 20) # hide
