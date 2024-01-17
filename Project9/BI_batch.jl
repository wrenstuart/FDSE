# Set orders of magnitude for Richardson numbers Ri and vertical shear parameters λ
log_Ri = [0]    # For Ri
log_α = [2]                # For s

include("BI.jl")

# Iterate over the orders of magnitude
for i₁ = 1 : length(log_Ri)
    Ri = 10 ^ log_Ri[i₁]        # Set Ri
    for i₂ = 1: length(log_α)
        α = 10 ^ log_α[i₂]     # Set λ
        label = "_Ri_" * string(log_Ri[i₁]) * "_α_" * string(log_α[i₂])
        @info label
        run_BI_sim(Ri, α, label)
    end
end