# This script simulates Kelvin-Helmholtz instability in 2D using Oceananigans

@info "starting"

# Load some standard libraries that we will need
using Printf
using Oceananigans
using Oceananigans.Units

# First, we need to set some physical parameters for the simulation
# Set the domain size in non-dimensional coordinates
Lx = 100kilometers  # size in the x-direction (east)
Ly = 300kilometers  # size in the y-direction (north)
Lz = 50meters       # size in the vertical (z) direction 

# Set the grid size
Nx = 32  # number of gridpoints in the x-direction
Ny = 128
Nz = 32   # number of gridpoints in the z-direction

# Some timestepping parameters
max_Δt = 1hour # maximum allowable timestep 
duration = 10days # The non-dimensional duration of the simulation

# Set the Reynolds number (Re=Ul/ν)
Re = 5000
# Set the Prandtl number (Pr=ν/κ)
Pr = 1

# Parameters for the initial condition:
ρ₀ = 1000   # reference density
f = 1e-4/second    # coriolis paramtere
S₀ = 25     # reference salinity
T₀ = 25     # reference temperature
N² = 1e-4  # vertical stratification due to temperature
M² = 1e-6  # horizontal stratification due to salinity
α = 0.17/ρ₀  # coefficient of thermal expansion
β = 0.76/ρ₀  # coefficient of haline contraction
g = 9.81



# Set the amplitude of the random perturbation (kick)
kick = 0.01

# construct a rectilinear grid using an inbuilt Oceananigans function
# Here, we use periodic (cyclic) boundary conditions in x
grid = RectilinearGrid(size = (Nx, Ny, Nz), x = (0, Lx), y = (0, Ly), z = (-Lz, 0), topology = (Periodic, Bounded, Bounded))

# No boundary conditions explicitly set - the BCs will default to free-slip and no-flux in z

# Now, define a 'model' where we specify the grid, advection scheme, bcs, and other settings
model = NonhydrostaticModel(; grid,
              advection = UpwindBiasedFifthOrder(),  # Specify the advection scheme.  Another good choice is WENO() which is more accurate but slower
            timestepper = :RungeKutta3, # Set the timestepping scheme, here 3rd order Runge-Kutta
                tracers = (:T, :S),  # Set the name(s) of any tracers, here S is salinity and T is temperature
               buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(thermal_expansion=α, haline_contraction=β)), # this tells the model that b will act as the buoyancy (and influence momentum) 
                closure = (ScalarDiffusivity(ν = 1 / Re, κ = 1 / Re)),  # set a constant kinematic viscosity and diffusivty, here just 1/Re since we are solving the non-dimensional equations 
                coriolis = FPlane(f=f)  # set coriolis parameter
)

# Set initial conditions
# start with linear front
#uᵢ(x, y, z) = -M²/f*(z+Lz) + kick*randn()
uᵢ(x, y, z) = M²/f*(z+Lz) + kick*randn()
vᵢ(x, y, z) = kick * randn()
wᵢ(x, y, z) = kick * randn()
Tᵢ(x, y, z) = T₀ + N²/(g*α)*z
Sᵢ(x, y, z) = S₀ + M²/(g*β)*y

# Send the initial conditions to the model to initialize the variables
set!(model, u = uᵢ, v = vᵢ, w = wᵢ, S = Sᵢ, T=Tᵢ)

# Now, we create a 'simulation' to run the model for a specified length of time
simulation = Simulation(model, Δt = max_Δt, stop_time = duration)

# ### The `TimeStepWizard`
wizard = TimeStepWizard(cfl = 0.85, max_change = 1.1, max_Δt = max_Δt)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

# ### A progress messenger
start_time = time_ns()
progress(sim) = @printf("i: % 6d, sim time: % 10s, wall time: % 10s, Δt: % 10s, CFL: %.2e\n",
                        sim.model.clock.iteration,
                        prettytime(sim.model.clock.time),
                        prettytime(1e-9 * (time_ns() - start_time)),
                        sim.Δt,
                        AdvectiveCFL(sim.Δt)(sim.model))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

# ### Output
u, v, w = model.velocities # unpack velocity `Field`s
S = model.tracers.S # extract the salinity
T = model.tracers.T # extract the temperature

# Now, calculate secondary quantities
# Oceananigans has functions to calculate derivatives on the model grid

ζ = ∂x(v)-∂y(u)

# Set the name of the output file
filename = "BI"

simulation.output_writers[:xy_slices] =
    JLD2OutputWriter(model, (; u, v, w, T, S, ζ),
                          filename = filename * ".jld2",
                          indices = (:, :, 1),
                         schedule = TimeInterval(1hour),
                            overwrite_existing = true)

nothing # hide

# Now, run the simulation
run!(simulation)

# After the simulation is different, plot the results and save a movie
include("plot_BI-john.jl")
