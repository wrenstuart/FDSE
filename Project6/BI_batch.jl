# An example of a 2D gravity current (lock-release problem) using Oceananigans

# Load some standard libraries that we will need
using Printf

M²s = [1.1e-6 1.5e-6 1e-6 9e-7 8e-7 7e-7 5e-7 3e-7 7e-7 1e-7 1e-8]
N²s = [1e-4 2.25e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 4.9e-5 1e-4 1e-4]
f = 1e-4

# giving M², N², f, simnum

for i in 1:length(M²s)
  if i > 10
    global M², N², simnum
    simnum = i
    M² = M²s[i]
    N² = N²s[i]
    include("BI.jl")
  end
end