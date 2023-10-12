# An example of a 2D gravity current (lock-release problem) using Oceananigans

# Load some standard libraries that we will need
using Printf
using Oceananigans
using Oceananigans.Units

# First, we need to set some physical parameters for the simulation
# Set the domain size in non-dimensional coordinates
Lx = 100kilometers  # size in the x-direction
Ly = 300kilometers  # size in the x-direction
Lz = 50   # size in the vertical (z) direction 

# Set the grid size
Nx = 64  # number of gridpoints in the x-direction
Ny = 256  # number of gridpoints in the x-direction
Nz = 32   # number of gridpoints in the z-direction

# Some timestepping parameters
max_Δt = 1hour # maximum allowable timestep 
duration = 10days # The non-dimensional duration of the simulation

# Set the amplitude of the random perturbation (kick)
kick = 0.01  # random velocty perturbation

# Now, some parameters that will be used for the initial conditions
xl = Lx / 10 # The location of the 'lock'
Lf = Lx / 100 # The width of the initial buoyancy step

f = 1e-4 # Coriolis parameter
g = 9.81 # gravitational acceleration
ρ₀ = 1020 # reference density

α = 0.17/ρ₀   #  Thermal expansion coefficient (TODO: check factor of 1000)
β = 0.76/ρ₀ # Haline contraction coefficient

M2 = 1e-6 # Initial horizontal buoyancy gradient
N2 = 1e-4 # Initial vertical buoyancy gradient

T₀ = 10 # Reference temperature (units = ᵒC) 
S₀ = 20 # Reference salinity (units = g/kg)

# construct a rectilinear grid using an inbuilt Oceananigans function
# Here, the topology parameter sets the style of boundaries in the x, y, and z directions
# 'Bounded' corresponds to wall-bounded directions and 'Flat' corresponds to the dimension that is not considered (here, that is the y direction)
grid = RectilinearGrid(size = (Nx, Ny, Nz), x = (0, Lx), y = (0, Ly), z = (-Lz, 0), topology = (Periodic, Bounded, Bounded))

# Now, define a 'model' where we specify the grid, advection scheme, bcs, and other settings
model = NonhydrostaticModel(; grid,
              advection = UpwindBiasedFifthOrder(),  # Specify the advection scheme.  Another good choice is WENO() which is more accurate but slower
            timestepper = :RungeKutta3, # Set the timestepping scheme, here 3rd order Runge-Kutta
                tracers = (:T, :S),  # Set the name(s) of any tracers, here b is buoyancy and c is a passive tracer (e.g. dye)
               buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(thermal_expansion=α, haline_contraction=β)), # this tells the model that b will act as the buoyancy (and influence momentum) 
               coriolis = coriolis = FPlane(f=f), # this line tells the mdoel not to include system rotation (no Coriolis acceleration)
)

# Set initial conditions
# Here, we start with a tanh function for buoyancy and add a random perturbation to the velocity. 
uᵢ(x, y, z) = M2/f*(z+Lz) + kick * randn()
vᵢ(x, y, z) = kick * randn() 
wᵢ(x, y, z) = kick * randn()
Tᵢ(x, y, z) = T₀ + N2/(g*α)*z
Sᵢ(x, y, z) = S₀ + M2/(g*β)*y

# Send the initial conditions to the model to initialize the variables
set!(model, u = uᵢ, v = vᵢ, w = wᵢ, T = Tᵢ, S = Sᵢ)

# Now, we create a 'simulation' to run the model for a specified length of time
simulation = Simulation(model, Δt = 1minute, stop_time = duration)

# ### The `TimeStepWizard`
#
# The TimeStepWizard manages the time-step adaptively, keeping the
# Courant-Freidrichs-Lewy (CFL) number close to `1.0` while ensuring
# the time-step does not increase beyond the maximum allowable value
wizard = TimeStepWizard(cfl = 0.85, max_change = 1.1, max_Δt = max_Δt)
# A "Callback" pauses the simulation after a specified number of timesteps and calls a function (here the timestep wizard to update the timestep)
# To update the timestep more or less often, change IterationInterval in the next line
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# ### A progress messenger
# We add a callback that prints out a helpful progress message while the simulation runs.

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
T = model.tracers.T # extract the buoyancy
S = model.tracers.S # extract the tracer

ζ = ∂x(v) - ∂y(u) # The vertical vorticity

# Set the name of the output file
filename = "BI_xz"

simulation.output_writers[:xz_slices] =
    JLD2OutputWriter(model, (; u, v, w, T, S, ζ),
                          filename = filename * ".jld2",
                          indices = (:, 1, :),
                         schedule = TimeInterval(1hour),
                            overwrite_existing = true)

# If you are running in 3D, you could save an xy slice like this:                             
filename = "BI_xy"
simulation.output_writers[:xy_slices] =
    JLD2OutputWriter(model, (; u, v, w, T, S, ζ),
                          filename = filename * ".jld2",
                          indices = (:,:,30),
                        schedule = TimeInterval(1hour),
                            overwrite_existing = true)

nothing # hide

# Now, run the simulation
run!(simulation)

# After the simulation is different, plot the results and save a movie
include("plot_BI.jl")
