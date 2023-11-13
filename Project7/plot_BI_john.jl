
# This script reads in output from gravitycurrent.jl, makes a plot, and saves an animation

using Oceananigans, JLD2, Plots, Printf
using Oceananigans.Units

# Define some constants (these should match BI.jl)
f = 1

# Set the filename (without the extension)
filename_xy = "Project7/raw-output/BI_xy"
filename_xz = "Project7/raw-output/BI_xz"
filename_yz = "Project7/raw-output/BI_yz"

# Read in the first iteration.  We do this to load the grid
# filename * ".jld2" concatenates the extension to the end of the filename
u_ic = FieldTimeSeries(filename_xy * ".jld2", "u", iterations = 0)
v_ic = FieldTimeSeries(filename_xy * ".jld2", "v", iterations = 0)
w_ic = FieldTimeSeries(filename_xy * ".jld2", "w", iterations = 0)
b_ic = FieldTimeSeries(filename_xy * ".jld2", "b", iterations = 0)
ζ_ic = FieldTimeSeries(filename_xy * ".jld2", "ζ", iterations = 0)

## Load in coordinate arrays
## We do this separately for each variable since Oceananigans uses a staggered grid
xu, yu, zu = nodes(u_ic)
xv, yv, zv = nodes(v_ic)
xw, yw, zw = nodes(w_ic)
xb, yb, zb = nodes(b_ic)
xζ, yζ, zζ = nodes(ζ_ic)

## Now, open the file with our data
file_xy = jldopen(filename_xy * ".jld2")
file_xz = jldopen(filename_xz * ".jld2")
file_yz = jldopen(filename_yz * ".jld2")

## Extract a vector of iterations
iterations = parse.(Int, keys(file_xz["timeseries/t"]))

@info "Making an animation from saved data..."

t_save = zeros(length(iterations))

# Here, we loop over all iterations
anim = @animate for (i, iter) in enumerate(iterations)

    @info "Drawing frame $i from iteration $iter..."

    b_xz = file_xz["timeseries/b/$iter"][:, 1, :];
    ζ_xz = file_xz["timeseries/ζ/$iter"][:, 1, :];

    b_xy = file_xy["timeseries/b/$iter"][:, :, 1];
    ζ_xy = file_xy["timeseries/ζ/$iter"][:, :, 1];

    global b_yz = file_yz["timeseries/b/$iter"][1, :, :];
    global v_yz = file_yz["timeseries/v/$iter"][1, :, :];
    global ζ_yz = file_yz["timeseries/ζ/$iter"][1, :, :];

    t = file_xz["timeseries/t/$iter"];

    # Save some variables to plot at the end
    t_save[i] = t # save the time

        b_xy_plot = heatmap(xb/1kilometer, yb/1kilometer, b_xy'; color = :haline, xlabel = "x (km)", ylabel = "y (km)");  
        ζ_xy_plot = heatmap(xζ/1kilometer, yζ/1kilometer, ζ_xy'/f; color = :balance, xlabel = "x (km)", ylabel = "y (km)", clims=(-0.5,0.5));  

        #ζ_yz_plot = heatmap(yζ/1kilometer, zζ, ζ_yz'/f; color = :balance, xlabel = "y (km)", ylabel = "z (m)", clims=(-3,3));
        ζ_yz_plot = heatmap(yζ/1kilometer, zζ, ζ_yz'/f.*20, clims=(-21,21), color = :balance);
        v_yz_plot = heatmap(yv/1kilometer, zv, v_yz', color = :balance);
        b_yz_plot = heatmap(yb/1kilometer, zb, b_yz', color = :haline);
        contour!(yb/1kilometer, zb, b_yz')#, levels=LinRange(20,20.6,50), color = :black)

    b_title = @sprintf("b");
    ζ_title = @sprintf("ζ/f");

# Combine the sub-plots into a single figure
    plot(b_xy_plot, ζ_xy_plot, layout = (1, 2), size = (1402, 600),
    title = [b_title ζ_title])
    #plot(b_yz_plot, size=(802,200))

    iter == iterations[end] && close(file_xz)
end

close(file_xy)
close(file_yz)
# Save the animation to a file
mp4(anim, "Project7/videos/BI.mp4", fps = 20) # hide
