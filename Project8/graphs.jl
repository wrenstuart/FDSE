using Oceananigans, JLD2, Plots, Printf
using Oceananigans.Units

function ℬ_timeseries(label)
    filename_ℬ = "Project8/raw-output/BI_ℬ" * label
    file_ℬ = jldopen(filename_ℬ * ".jld2")
    ℬ_ic = FieldTimeSeries(filename_ℬ * ".jld2", "ℬ", iterations = 0)
    xℬ, yℬ, zℬ = nodes(ℬ_ic)
    iterations = parse.(Int, keys(file_ℬ["timeseries/t"]))
    t_save = zeros(length(iterations))
    ℬ_mean = zeros(Float64, length(iterations))
    for i = 11 : length(iterations)
        iter = iterations[i]
        t = file_ℬ["timeseries/t/$iter"]
        t_save[i] = t
        if i > 10       # Ignore the first few iterations due to noise
            ℬ_xy = file_ℬ["timeseries/ℬ/$iter"][:, 1, :]
            ℬ_mean[i] = sum(ℬ_xy) / Float64(length(xℬ) * length(yℬ))
        end
    end

    return t_save, ℬ_mean

end

function ζ_var_timeseries(label)
    filename_xy = "Project8/raw-output/BI_xy" * label
    file_xy = jldopen(filename_xy * ".jld2")
    ζ_ic = FieldTimeSeries(filename_xy * ".jld2", "ζ", iterations = 0)
    xζ, yζ, zζ = nodes(ζ_ic)
    iterations = parse.(Int, keys(file_xy["timeseries/t"]))
    t_save = zeros(length(iterations))
    ζ_var = zeros(Float64, length(iterations))
    for i = 11 : length(iterations)
        iter = iterations[i]
        t = file_xy["timeseries/t/$iter"]
        t_save[i] = t
        if i > 10       # Ignore the first few iterations due to noise
            ζ_xy = file_xy["timeseries/ζ/$iter"]
            ζ_var[i] = sum(ζ_xy .* ζ_xy) / Float64(length(xζ) * length(yζ))
        end
    end

    return t_save, ζ_var

end

function δ_var_timeseries(label)
    filename_xy = "Project8/raw-output/BI_xy" * label
    file_xy = jldopen(filename_xy * ".jld2")
    δ_ic = FieldTimeSeries(filename_xy * ".jld2", "δ", iterations = 0)
    xδ, yδ, zδ = nodes(δ_ic)
    iterations = parse.(Int, keys(file_xy["timeseries/t"]))
    t = zeros(length(iterations))
    δ_var = zeros(Float64, length(iterations))
    for i = 11 : length(iterations)
        iter = iterations[i]
        t[i] = file_xy["timeseries/t/$iter"]
        if i > 10       # Ignore the first few iterations due to noise
            δ_xy = file_xy["timeseries/δ/$iter"]
            δ_var[i] = sum(δ_xy .* δ_xy) / Float64(length(xδ) * length(yδ))
        end
    end

    return t, δ_var

end

function δ_skew_timeseries(label)
    filename_xy = "Project8/raw-output/BI_xy" * label
    file_xy = jldopen(filename_xy * ".jld2")
    δ_ic = FieldTimeSeries(filename_xy * ".jld2", "δ", iterations = 0)
    xδ, yδ, zδ = nodes(δ_ic)
    iterations = parse.(Int, keys(file_xy["timeseries/t"]))
    t = zeros(length(iterations))
    δ_skew = zeros(Float64, length(iterations))
    for i = 11 : length(iterations)
        iter = iterations[i]
        t[i] = file_xy["timeseries/t/$iter"]
        δ_xy = file_xy["timeseries/δ/$iter"]
        δ_var = sum(δ_xy .^ 2) / Float64(length(xδ) * length(yδ))
        δ_scale = δ_var ^ 0.5
        δ_skew[i] = sum(δ_xy .^ 3) / (length(xδ) * length(yδ) * δ_scale ^ 3)
    end

    return t, δ_skew

end

function ζ_skew_timeseries(label)
    filename_xy = "Project8/raw-output/BI_xy" * label
    file_xy = jldopen(filename_xy * ".jld2")
    ζ_ic = FieldTimeSeries(filename_xy * ".jld2", "ζ", iterations = 0)
    xζ, yζ, zζ = nodes(ζ_ic)
    iterations = parse.(Int, keys(file_xy["timeseries/t"]))
    t_save = zeros(length(iterations))
    ζ_skew= zeros(Float64, length(iterations))
    for i = 11 : length(iterations)
        iter = iterations[i]
        t = file_xy["timeseries/t/$iter"]
        t_save[i] = t
        ζ_xy = file_xy["timeseries/ζ/$iter"]
        ζ_var = sum(ζ_xy .^ 2) / (length(xζ) * length(yζ))
        ζ_scale = ζ_var ^ 0.5
        ζ_skew[i] = sum(ζ_xy .^ 3) / (length(xζ) * length(yζ) * ζ_scale ^ 3)
    end

    return t_save, ζ_skew

end

#log_Ri = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]    # For Ri
log_Ri = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
log_α = [2]                # For s

include("BI.jl")

f = 1e-4
#i_break = 0
for i₁ = 1 : length(log_Ri)
    for i₂ = 1: length(log_α)
        #=if isinteger(log_Ri[i₁])
            Ri_label = Int64(log_Ri[i₁])
        else
            Ri_label = log_Ri[i₁]
        end
        =#
        Ri_label = log_Ri[i₁]
        label = "_Ri_" * string(Ri_label) * "_α_" * string(log_α[i₂])
        t, ζ_var = ζ_var_timeseries(label)
        ~, δ_var = δ_var_timeseries(label)
        ~, ℬ_mean = ℬ_timeseries(label)
        ~, ζ_skew = ζ_skew_timeseries(label)
        ~, δ_skew = δ_skew_timeseries(label)
        t = t/t[end] * 40
        plt = Plots.plot(t, δ_var/maximum(δ_var))
        Plots.plot!(t, ζ_var/maximum(ζ_var))
        Plots.plot!(t, ℬ_mean/maximum(ℬ_mean))
        Plots.plot!(t, ζ_skew)
        Plots.plot!(t, δ_skew/maximum(-δ_skew))
        #Plots.plot!(t, δ_skew)
        display(plt)
        #=global i_break = length(t)
        for i = Int64(round(length(t)/2)) : length(t)
            if ζ_skew[i] < ζ_skew[i-1]
                global i_break = i-1
                @info i_break
                break
            end
        end=#

        ζ_max_skew = maximum(ζ_skew)
        if Ri_label > 2
            ζ_max_skew = maximum(ζ_skew[Int64(round(length(t) * 0.5)) : Int64(round(length(t) * 0.75))])
        end
        i_break = findall(x -> x == ζ_max_skew, ζ_skew)[1]
        if Ri_label > 3
            δ_max_skew = maximum( - δ_skew[Int64(round(length(t) * 0.5)) : Int64(round(length(t) * 0.55))])
            i_break = findall(x -> x == -δ_max_skew, δ_skew)[1]
            ζ_max_skew = ζ_skew[i_break]
        end
        
        #@info i_break, t[i_break]
        filename_xy = "Project8/raw-output/BI_xy" * label
        file_xy = jldopen(filename_xy * ".jld2")
        ζ_ic = FieldTimeSeries(filename_xy * ".jld2", "ζ", iterations = 0)
        δ_ic = FieldTimeSeries(filename_xy * ".jld2", "δ", iterations = 0)
        iterations = parse.(Int, keys(file_xy["timeseries/t"]))
        
        iter = iterations[i_break]
        xζ, yζ, zζ = nodes(ζ_ic)
        xδ, yδ, zδ = nodes(δ_ic)
        δ_xy = file_xy["timeseries/δ/$iter"][:, :, 1];
        ζ_xy = file_xy["timeseries/ζ/$iter"][:, :, 1];
        mean_δ = sum(δ_xy) / (length(xδ) * length(yδ))
        mean_ζ = sum(ζ_xy) / (length(xζ) * length(yζ))
        δ_scale = δ_var[i_break] ^ 0.5
        ζ_scale = ζ_var[i_break] ^ 0.5
        mean_δ³ = sum((δ_xy .- mean_δ) .^ 3) / (length(xδ) * length(yδ))
        mean_ζ³ = sum((ζ_xy .- mean_ζ) .^ 3) / (length(xζ) * length(yζ))
        skew_δ = mean_δ³ / (δ_scale ^ 3)
        skew_ζ = mean_ζ³ / (ζ_scale ^ 3)
        @info Ri_label, skew_ζ, skew_δ
        plt = Plots.scatter(ζ_xy/f, δ_xy/f, primary = false, markeralpha = 0.01, markerstrokewidth = 0, color = :black, xlims = (-5*ζ_scale/f, 5*ζ_scale/f), ylims = (-5*δ_scale/f, 5*δ_scale/f))
        #plt = scatter(ζ_xy/f, δ_xy/f, primary = false, markeralpha = 0.02, markerstrokewidth = 0, color = :black, xlims = (-5, 5), ylims = (-5, 5), framestyle = :origin)
        Plots.plot!([-1, -1], [-5, 5], linestyle = :dash, color = :red, legend = false)
        xlabel!("\$ζ/f\$")
        ylabel!("\$δ/f\$")
        title!("\$\\mathrm{Ri}=" * string(round(10^Ri_label, sigdigits = 3)) * "\$")
        annotate!(3.5*ζ_scale*1e4, 3.5*δ_scale*1e4, "\$\\tilde\\mu_{3,\\zeta}=" * string(round(10^skew_δ, sigdigits = 3)) * "\$\n\$\\tilde\\mu_{3,\\delta}=" * string(round(10^skew_δ, sigdigits = 3)) * "\$")
        display(plt)
        
    end
end

