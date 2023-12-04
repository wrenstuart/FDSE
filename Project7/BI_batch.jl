# Set orders of magnitude for Richardson numbers Ri and vertical shear parameters λ
S₁ = [0, 1, 2, 3, 4]             # For Ri
S₂ = [1]          # For λ

include("BI.jl")

# Iterate over the orders of magnitude
for i₁ = 1 : length(S₁)
    Ri = 10 ^ S₁[i₁]        # Set Ri
    for i₂ = 1: length(S₂)
        λ = 10 ^ S₂[i₂]     # Set λ
        label = "_Ri-" * string(S₁[i₁]) * "_λ-" * string(S₂[i₂])
        @info label
        run_BI_sim(Ri, λ, label)
    end
end