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

Ro = log_Ri * 0;
δ_on_f = log_Ri * 0;
ζ_skew_of_Ri = log_Ri * 0;
δ_skew_of_Ri = log_Ri * 0;

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
        display(plt)

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
        iterations = parse.(Int, keys(file_xy["timeseries/t"]))
        
        iter = iterations[i_break]
        δ_xy = file_xy["timeseries/δ/$iter"][:, :, 1];
        ζ_xy = file_xy["timeseries/ζ/$iter"][:, :, 1];
        δ_scale = δ_var[i_break] ^ 0.5
        ζ_scale = ζ_var[i_break] ^ 0.5
        @info Ri_label, ζ_skew[i_break], δ_skew[i_break]
        plt = Plots.scatter(ζ_xy/f, δ_xy/f, primary = false, markeralpha = 0.01, markerstrokewidth = 0, color = :black, xlims = (-5*ζ_scale/f, 5*ζ_scale/f), ylims = (-5*δ_scale/f, 5*δ_scale/f))
        Plots.plot!([-1, -1], [-5, 5], linestyle = :dash, color = :red, legend = false)
        Plots.xlabel!("\$ζ/f\$")
        Plots.ylabel!("\$δ/f\$")
        title!("\$\\mathrm{Ri}=" * string(round(10^Ri_label, sigdigits = 3)) * "\$")
        annotate!(2*ζ_scale*1e4, 2*δ_scale*1e4, "\$\\tilde\\mu_{3,\\zeta}=" * string(round(ζ_skew[i_break], sigdigits = 3)) * "\$\n\$\\tilde\\mu_{3,\\delta}=" * string(round(δ_skew[i_break], sigdigits = 3)) * "\$", :left)
        display(plt)
        savefig("Project8/graphs/snapshot" * label * ".png")
        
        Ro[i₁] = ζ_scale * 1e4
        δ_on_f[i₁] = δ_scale * 1e4
        ζ_skew_of_Ri[i₁] = ζ_skew[i_break]
        δ_skew_of_Ri[i₁] = δ_skew[i_break]

    end
end

Ri = 10 .^ log_Ri
plt = Plots.plot(log10.(1 .+ Ri), log10.(Ro))
Plots.plot!([0, 4], [0.1, -0.9], linestyle = :dash, color = :black, legend = false)
Plots.xlabel!("\$\\log(1+\\mathrm{Ri})\$")
Plots.ylabel!("\$\\log(\\mathrm{Ro})\$")
display(plt)
savefig("Project8/graphs/ζ_Ri.png")
plt = Plots.plot(log10.(1 .+ Ri), log10.(δ_on_f))
Plots.plot!([0, 4], [-0.2, -2.2], linestyle = :dash, color = :black, legend = false)
Plots.xlabel!("\$\\log(1+\\mathrm{Ri})\$")
Plots.ylabel!("\$\\frac{1}{2}\\log(\\langle\\delta^2\\rangle/f^2)\$")
display(plt)
savefig("Project8/graphs/δ_Ri.png")
plt = Plots.plot(log10.(1 .+ Ri), log10.(Ro ./ δ_on_f), aspect_ratio = :equal)
Plots.plot!([0, 4], [0.3, 1.3], linestyle = :dash, color = :black, legend = false)
Plots.xlabel!("\$\\log(1+\\mathrm{Ri})\$")
Plots.ylabel!("\$\\frac{1}{2}\\log(\\langle\\zeta^2\\rangle/\\langle\\delta^2\\rangle)\$")
display(plt)
savefig("Project8/graphs/ζ-on-δ_Ri.png")
plt = Plots.plot(log10.(1 .+ Ri), ζ_skew_of_Ri, legend = false)
Plots.xlabel!("\$\\log(1+\\mathrm{Ri})\$")
Plots.ylabel!("\$\\tilde\\mu_{3,\\zeta}\$")
display(plt)
savefig("Project8/graphs/ζ-skew_Ri.png")
plt = Plots.plot(log10.(1 .+ Ri), δ_skew_of_Ri, legend = false)
Plots.xlabel!("\$\\log(1+\\mathrm{Ri})\$")
Plots.ylabel!("\$\\tilde\\mu_{3,\\delta}\$")
display(plt)