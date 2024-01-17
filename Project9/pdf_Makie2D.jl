# This script reads in output from KH.jl, makes a plot, and saves an animation

using Oceananigans, JLD2, Printf, CairoMakie

CairoMakie.activate!()

pre = "Project9/raw-output/BI_xy"
Ris = [0, 1, 2, 3, 4]
post = ""

n_nodes = 100
freq = zeros(n_nodes, n_nodes)
fixed_axes = false

function ζ_var_timeseries(label)
    filename_xy = "Project9/raw-output/BI_xy" * label
    file_xy = jldopen(filename_xy * ".jld2")
    iterations = parse.(Int, keys(file_xy["timeseries/t"]))
    t_save = zeros(length(iterations))
    ζ_ic = FieldTimeSeries(filename_xy * ".jld2", "ζ", iterations = 0)
    xζ, yζ, zζ = nodes(ζ_ic)
    ζ_var_mean = zeros(Float64, length(iterations))
    for i = 11 : length(iterations)
        iter = iterations[i]
        t = file_xy["timeseries/t/$iter"]
        t_save[i] = t
        if i > 10       # Ignore the first few iterations due to noise
            ζ_xy = file_xy["timeseries/ζ/$iter"]
            ζ_var_mean[i] = sum(ζ_xy .* ζ_xy) / Float64(length(xζ) * length(yζ))
        end
    end

    return t_save, ζ_var_mean

end

function δ_var_timeseries(label)
    filename_xy = "Project9/raw-output/BI_xy" * label
    file_xy = jldopen(filename_xy * ".jld2")
    iterations = parse.(Int, keys(file_xy["timeseries/t"]))
    t = zeros(length(iterations))
    δ_ic = FieldTimeSeries(filename_xy * ".jld2", "δ", iterations = 0)
    xδ, yδ, zδ = nodes(δ_ic)
    δ_var_mean = zeros(Float64, length(iterations))
    for i = 11 : length(iterations)
        iter = iterations[i]
        t[i] = file_xy["timeseries/t/$iter"]
        if i > 10       # Ignore the first few iterations due to noise
            δ_xy = file_xy["timeseries/δ/$iter"]
            δ_var_mean[i] = sum(δ_xy .* δ_xy) / Float64(length(xδ) * length(yδ))
        end
    end

    return t, δ_var_mean

end

for i = 1 : length(Ris)

    log_Ri = Ris[i]
    Ri = 10 ^ log_Ri

    log_α = 2

    #=if isinteger(log_Ri)
        Ri_label = Int64(log_Ri)
    else
        Ri_label = log_Ri
    end
    =#
    Ri_label = log_Ri

    label = "_Ri_" * string(Ri_label) * "_α_" * string(log_α)

    # Set the filename (without the extension)
    filename = pre * label

    # Read in the first iteration.  We do this to load the grid
    # filename * ".jld2" concatenates the extension to the end of the filename
    ζ_ic = FieldTimeSeries(filename * ".jld2", "ζ", iterations = 0)
    δ_ic = FieldTimeSeries(filename * ".jld2", "δ", iterations = 0)

    ## Load in coordinate arrays
    ## We do this separately for each variable since Oceananigans uses a staggered grid
    xζ, yζ, zζ = nodes(ζ_ic)
    xδ, yδ, zδ = nodes(δ_ic)

    ## Now, open the file with our data
    file_xy = jldopen(filename * ".jld2")

    ## Extract a vector of iterations
    iterations = parse.(Int, keys(file_xy["timeseries/t"]))

    @info "Making an animation from saved data..."
    
    t, ζ_var_mean = ζ_var_timeseries(label)
    ~, δ_var_mean = δ_var_timeseries(label)

    if fixed_axes
        ζ_lim = 5
        δ_lim = 5
    else
        ζ_lim = 5 * maximum(ζ_var_mean) ^ 0.5 * 1e4
        δ_lim = 5 * maximum(δ_var_mean) ^ 0.5 * 1e4
    end

    f = Figure(backgroundcolor=:white, reolsution=(800,300))
    ax = Axis(f[1,1], limits = ((-ζ_lim, ζ_lim), (-δ_lim, δ_lim)))

    function update_frame(iter)
        i, time = iter;
        i = i + 2
        empty!(ax)
        global freq = zeros(n_nodes, n_nodes)

        x = collect(range(-ζ_lim + ζ_lim/(2*n_nodes), ζ_lim - ζ_lim/(2*n_nodes), n_nodes))
        y = collect(range(-δ_lim + δ_lim/(2*n_nodes), δ_lim - δ_lim/(2*n_nodes), n_nodes))

        for j = i - 2 : i + 2
        # CHANGE LINE BELOW IF CHANGE LINE ABOVE
        time = iterations[j+Int64(round(length(iterations)/2))]
        ζ_xy = file_xy["timeseries/ζ/$time"];
        δ_xy = file_xy["timeseries/δ/$time"];
        diff_ζ = 2 * ζ_lim / n_nodes
        diff_δ = 2 * δ_lim / n_nodes

        for row = 1 : size(ζ_xy)[1]
            for col = 1 : size(ζ_xy)[2]
                ζ = ζ_xy[row, col]
                δ = δ_xy[row, col]
                x_box = Int64(ceil((ζ_lim + ζ*1e4) / diff_ζ))
                y_box = Int64(ceil((δ_lim + δ*1e4) / diff_δ))
                if 0 < x_box && x_box <= n_nodes && 0 < y_box && y_box <= n_nodes
                    global freq[x_box, y_box] += 1
                end
            end
        end
        
        end

        CairoMakie.heatmap!(x, y, -log.(1 .+ freq), colormap = :bilbao)
        vlines!(ax, 0, color = :black)
        hlines!(ax, 0, color = :black)
        hidespines!(ax)

    end

    if fixed_axes
        fix = "fixed-axes"
    else
        fix = "variable-axes"
    end

    record(update_frame, f, "Project9/histograms/hist_" * fix * label * ".mp4", enumerate(iterations[3+Int64(round(length(iterations)/2)):1:end-3]); framerate = 20)

    close(file_xy)
    nothing

end