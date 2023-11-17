# Eady model

# Load some standard libraries that we will need
using Printf
using Oceananigans
using Oceananigans.Units

# Set the grid size
Nx = 128  # number of gridpoints in the x-direction
Ny = 128  # number of gridpoints in the x-direction
Nz = 16   # number of gridpoints in the z-direction

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

# Set some timestepping parameters
max_Δt = T / 10
duration = T * 40
@printf("Simulation will last %s", prettytime(duration))

# Set relative amplitude for random velocity perturbation
kick = 0.05

# Define the background fields
B₀(x, y, z, t) = M² * y + N² * z   # Buoyancy
U₀(x, y, z, t) = -M²/f * (z + Lz)  # Velocity (initialised in thermal wind balance)
B_field = BackgroundField(B₀)
U_field = BackgroundField(U₀)

# construct a rectilinear grid using an inbuilt Oceananigans function
# Here, the topology parameter sets the style of boundaries in the x, y, and z directions
grid = RectilinearGrid(size = (Nx, Ny, Nz), x = (0, Lx), y = (0, Ly), z = (-Lz, 0), topology = (Periodic, Periodic, Bounded))

# Leave Oceananigans to set default boundary conditions for the perturbation fields

# Now, define a 'model' where we specify the grid, advection scheme, bcs, and other settings
model = NonhydrostaticModel(; grid,
              advection = UpwindBiasedFifthOrder(),  # Specify the advection scheme.  Another good choice is WENO() which is more accurate but slower
              timestepper = :RungeKutta3, # Set the timestepping scheme, here 3rd order Runge-Kutta
              tracers = :b,  # Set the name(s) of any tracers, here b is buoyancy and c is a passive tracer (e.g. dye)
              buoyancy = Buoyancy(model=BuoyancyTracer()), # this tells the model that b will act as the buoyancy (and influence momentum)
              background_fields = (b = B_field, u = U_field),
              coriolis = coriolis = FPlane(f = f)
)

# Set initial conditions, a random velocity perturbation
uᵢ(x, y, z) = U * kick * randn()
vᵢ(x, y, z) = U * kick * randn() 
wᵢ(x, y, z) = U * kick * randn()
bᵢ(x, y, z) = 0

# Send the initial conditions to the model to initialize the variables
set!(model, u = uᵢ, v = vᵢ, w = wᵢ, b = bᵢ)

# Now, we create a 'simulation' to run the model for a specified length of time
simulation = Simulation(model, Δt = max_Δt/10, stop_time = duration)

# ### The `TimeStepWizard`
#
# The TimeStepWizard manages the time-step adaptively, keeping the
# Courant-Freidrichs-Lewy (CFL) number close to `1.0` while ensuring
# the time-step does not increase beyond the maximum allowable value
wizard = TimeStepWizard(cfl = 0.85, max_change = 1.1, max_Δt = max_Δt)
# A "Callback" pauses the simulation after a specified number of timesteps and calls a function (here the timestep wizard to update the timestep)
# To update the timestep more or less often, change IterationInterval in the next line
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

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

u = Field(model.velocities.u + model.background_fields.velocities.u) # unpack velocity `Field`s
v = model.velocities.v
w = model.velocities.w
b = Field(model.tracers.b + model.background_fields.tracers.b) # extract the buoyancy and add the background field

ζ = ∂x(v) - ∂y(u)   # The vertical vorticity
δ = ∂x(u) + ∂y(v)   # The horizontal divergence

# Output the slice y =0
filename = "Project7/raw-output/BI_xz"
simulation.output_writers[:xz_slices] =
    JLD2OutputWriter(model, (; u, v, w, b, ζ, δ),
                            filename = filename * ".jld2",
                            indices = (:, 1, :),
                            schedule = TimeInterval(T/10),
                            overwrite_existing = true)

# Output the slice z = 0
filename = "Project7/raw-output/BI_xy"
simulation.output_writers[:xy_slices] =
    JLD2OutputWriter(model, (; u, v, w, b, ζ, δ),
                            filename = filename * ".jld2",
                            indices = (:, :, Nz),
                            schedule = TimeInterval(T/10),
                            overwrite_existing = true)

# Output the slice x = 0
filename = "Project7/raw-output/BI_yz"
simulation.output_writers[:yz_slices] =
    JLD2OutputWriter(model, (; u, v, w, b, ζ, δ),
                            filename = filename * ".jld2",
                            indices = (1, :, :),
                            schedule = TimeInterval(T/10),
                            overwrite_existing = true)

nothing # hide

# Now, run the simulation
run!(simulation)

# After the simulation is done, plot the results and save a movie
include("plot_BI_john.jl")