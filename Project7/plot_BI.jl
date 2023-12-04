using Oceananigans, JLD2, Plots, Printf
using Oceananigans.Units

# This script reads in output from gravitycurrent.jl, makes a plot, and saves an animation
function BI_plot(label)
# Set the two non-dimensional parameters
Ri = 100  # Richardson number, N²f²/M⁴
α = 10    # M²/f²

# Set the two dimensional parameters
H = 50    # Depth of mixed layer
f = 1e-4  # Coriolis parameter

# Calculate the other non-dimensional parameters
β = Ri * α^2    # N²/f²
Ro = Ri^(-0.5)  # (initial) Rossby number

# Calculate the other dimensional parameters and scales
M² = α * f^2                # Horizontal buoyancy gradient
N² = β * f^2                # Vertical buoyancy gradient
L = α * H * (1+Ri)^0.5      # Horizontal lengthscale (Rossby deformation radius)
U = α * f * H               # Velocity scale
T = (1+Ri)^0.5 / f          # Timescale of growth (and expected timescale of later flow, if Ro ∼ Ri^(-0.5)

# Set the domain size
Lx = 4 * 2*pi * L * 0.4^0.5 # Zonal extent, set to 4 wavelengths of the most unstable mode
Ly = (1+5^0.5)/2 * Lx       # Meridional extent, chosen to be an irrational multiple of Lx to avoid resonance
Lz = H                      # Vertical extent



# Set the filename (without the extension)
filename_xy = "Project7/raw-output/BI_xy" * label
filename_xz = "Project7/raw-output/BI_xz" * label
filename_yz = "Project7/raw-output/BI_yz" * label

# Read in the first iteration.  We do this to load the grid
# filename * ".jld2" concatenates the extension to the end of the filename
u_ic = FieldTimeSeries(filename_xy * ".jld2", "u", iterations = 0)
v_ic = FieldTimeSeries(filename_xy * ".jld2", "v", iterations = 0)
w_ic = FieldTimeSeries(filename_xy * ".jld2", "w", iterations = 0)
b_ic = FieldTimeSeries(filename_xy * ".jld2", "b", iterations = 0)
ζ_ic = FieldTimeSeries(filename_xy * ".jld2", "ζ", iterations = 0)
δ_ic = FieldTimeSeries(filename_xy * ".jld2", "δ", iterations = 0)


## Load in coordinate arrays
## We do this separately for each variable since Oceananigans uses a staggered grid
xu, yu, zu = nodes(u_ic)
xv, yv, zv = nodes(v_ic)
xw, yw, zw = nodes(w_ic)
xb, yb, zb = nodes(b_ic)
xζ, yζ, zζ = nodes(ζ_ic)
xδ, yδ, zδ = nodes(δ_ic)

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

    δ_xz = file_xz["timeseries/δ/$iter"][:, 1, :];
    ζ_xz = file_xz["timeseries/ζ/$iter"][:, 1, :];

    δ_xy = file_xy["timeseries/δ/$iter"][:, :, 1];
    ζ_xy = file_xy["timeseries/ζ/$iter"][:, :, 1];

    global δ_yz = file_yz["timeseries/δ/$iter"][1, :, :];
    global v_yz = file_yz["timeseries/v/$iter"][1, :, :];
    global ζ_yz = file_yz["timeseries/ζ/$iter"][1, :, :];

    t = file_xz["timeseries/t/$iter"];

    # Save some variables to plot at the end
    t_save[i] = t # save the time

        δ_xy_plot = heatmap(xδ/1kilometer, yδ/1kilometer, δ_xy'/f; color = :balance, xlabel = "x (km)", ylabel = "y (km)", clims=(-0.5,0.5));  
        ζ_xy_plot = heatmap(xζ/1kilometer, yζ/1kilometer, ζ_xy'/f; color = :balance, xlabel = "x (km)", ylabel = "y (km)", clims=(-0.5,0.5));  

        #ζ_yz_plot = heatmap(yζ/1kilometer, zζ, ζ_yz'/f; color = :balance, xlabel = "y (km)", ylabel = "z (m)", clims=(-3,3));
        ζ_yz_plot = heatmap(yζ/1kilometer, zζ, ζ_yz'/f.*20, clims=(-21,21), color = :balance);
        v_yz_plot = heatmap(yv/1kilometer, zv, v_yz', color = :balance);
        δ_yz_plot = heatmap(yδ/1kilometer, zδ, δ_yz', color = :haline);
        contour!(yδ/1kilometer, zδ, δ_yz')#, levels=LinRange(20,20.6,50), color = :black)

    δ_title = @sprintf("δ/f");
    ζ_title = @sprintf("ζ/f");

# Combine the sub-plots into a single figure
    plot(δ_xy_plot, ζ_xy_plot, layout = (1, 2), size = (1402, 600),
    title = [δ_title ζ_title])
    #plot(b_yz_plot, size=(802,200))

    iter == iterations[end] && close(file_xz)
end

close(file_xy)
close(file_yz)
# Save the animation to a file
mp4(anim, "Project7/videos/BI" * label * ".mp4", fps = 20) # hide
end