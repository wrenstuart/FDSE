using Oceananigans, JLD2, Plots, Printf
using Oceananigans.Units



function BI_plot(label)

# Set the two dimensional parameters
H = 50    # Depth of mixed layer
f = 1e-4  # Coriolis parameter

# Set the filename (without the extension)
filename_xy = "Project8/raw-output/BI_xy" * label
filename_xz = "Project8/raw-output/BI_xz" * label
filename_yz = "Project8/raw-output/BI_yz" * label
filename_ℬ = "Project8/raw-output/BI_ℬ" * label

# Read in the first iteration.  We do this to load the grid
u_ic = FieldTimeSeries(filename_xy * ".jld2", "u", iterations = 0)
v_ic = FieldTimeSeries(filename_xy * ".jld2", "v", iterations = 0)
w_ic = FieldTimeSeries(filename_xy * ".jld2", "w", iterations = 0)
b_ic = FieldTimeSeries(filename_xy * ".jld2", "b", iterations = 0)
ζ_ic = FieldTimeSeries(filename_xy * ".jld2", "ζ", iterations = 0)
δ_ic = FieldTimeSeries(filename_xy * ".jld2", "δ", iterations = 0)
ℬ_ic = FieldTimeSeries(filename_ℬ * ".jld2", "ℬ", iterations = 0)

# Load in co-ordinate arrays
# We do this separately for each variable since Oceananigans uses a staggered grid
xu, yu, zu = nodes(u_ic)
xv, yv, zv = nodes(v_ic)
xw, yw, zw = nodes(w_ic)
xb, yb, zb = nodes(b_ic)
xζ, yζ, zζ = nodes(ζ_ic)
xδ, yδ, zδ = nodes(δ_ic)
xℬ, yℬ, zℬ = nodes(ℬ_ic)

# Now, open the file with our data
file_xy = jldopen(filename_xy * ".jld2")
file_xz = jldopen(filename_xz * ".jld2")
file_yz = jldopen(filename_yz * ".jld2")
file_ℬ = jldopen(filename_ℬ * ".jld2")

# Extract a vector of iterations
iterations = parse.(Int, keys(file_xz["timeseries/t"]))
# iterations[i] is the number of timesteps passed for the i-th frame

@info "Making an animation from saved data..."

t_save = zeros(length(iterations))  # Contains the actual time

# Calculate the maximum relative vorticity and buoyancy flux to set the scale for the colourmap
ζ_max = 0
ℬ_max = 0
ℬ_mean = zeros(Float64, length(iterations))

for i = 11 : length(iterations)
    iter = iterations[i]
    t = file_xy["timeseries/t/$iter"]
    t_save[i] = t
    if i > 10       # Ignore the first few iterations due to noise
        ζ_xy = file_xy["timeseries/ζ/$iter"][:, 1, :]
        ℬ_xy = file_ℬ["timeseries/ℬ/$iter"][:, 1, :]
        this_ζ_max = maximum(ζ_xy)
        this_ℬ_max = maximum(ℬ_xy)
        ζ_max = maximum([ζ_max, this_ζ_max])
        ℬ_max = maximum([ℬ_max, this_ℬ_max])
        ℬ_mean[i] = sum(ℬ_xy) / Float64(length(xℬ) * length(yℬ))
    end
end

ℬ_mean_max, frame_ℬ_mean_max = findmax(ℬ_mean)
@info frame_ℬ_mean_max, frame_ℬ_mean_max/800

# Here, we loop over all iterations
anim = @animate for (i, iter) in enumerate(iterations)

    #@info "Drawing frame $i from iteration $iter..."

    #δ_xz = file_xz["timeseries/δ/$iter"][:, 1, :];
    #ζ_xz = file_xz["timeseries/ζ/$iter"][:, 1, :];

    δ_xy = file_xy["timeseries/δ/$iter"][:, :, 1];
    ζ_xy = file_xy["timeseries/ζ/$iter"][:, :, 1];
    ℬ_xy = file_ℬ["timeseries/ℬ/$iter"][:, :, 1];

    t = file_xy["timeseries/t/$iter"];  # The time that has elapsed by this iteration

    # Save some variables to plot at the end
    t_save[i] = t # save the time

        δ_xy_plot = Plots.heatmap(xδ/1kilometer, yδ/1kilometer, δ_xy'/f; color = :balance, xlabel = "x (km)", ylabel = "y (km)", clims=(-ζ_max/f, ζ_max/f));  
        ζ_xy_plot = Plots.heatmap(xζ/1kilometer, yζ/1kilometer, ζ_xy'/f; color = :balance, xlabel = "x (km)", ylabel = "y (km)", clims=(-ζ_max/f, ζ_max/f));  
        ℬ_xy_plot = Plots.heatmap(xℬ/1kilometer, yℬ/1kilometer, ℬ_xy'; color = :balance, xlabel = "x (km)", ylabel = "y (km)", clims=(-ℬ_max, ℬ_max));  

        #ζ_yz_plot = heatmap(yζ/1kilometer, zζ, ζ_yz'/f; color = :balance, xlabel = "y (km)", ylabel = "z (m)", clims=(-3,3));
        #ζ_yz_plot = heatmap(yζ/1kilometer, zζ, ζ_yz'/f.*20, clims=(-21,21), color = :balance);
        #v_yz_plot = heatmap(yv/1kilometer, zv, v_yz', color = :balance);
        #δ_yz_plot = heatmap(yδ/1kilometer, zδ, δ_yz', color = :haline);
        #contour!(yδ/1kilometer, zδ, δ_yz')#, levels=LinRange(20,20.6,50), color = :black)

    δ_title = @sprintf("δ/f");
    ζ_title = @sprintf("ζ/f");
    ℬ_title = @sprintf("Buoyancy flux ℬ")

    # Combine the sub-plots into a single figure
    Plots.plot(δ_xy_plot, ζ_xy_plot, ℬ_xy_plot, layout = (1, 3),# size = (1402, 600),
    title = [δ_title ζ_title ℬ_title])
    #plot(b_yz_plot, size=(802,200))

    iter == iterations[end] && close(file_xz)
end

close(file_xy)
close(file_yz)

# Save the animation to a file
mp4(anim, "Project8/videos/BI" * label * ".mp4", fps = 20) # hide

end