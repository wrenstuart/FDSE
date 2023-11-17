# Eady model

# Load some standard libraries that we will need
using Printf
using Oceananigans
using Oceananigans.Units


# Set the grid size
Nx = 128  # number of gridpoints in the x-direction
Ny = 128  # number of gridpoints in the x-direction
Nz = 16   # number of gridpoints in the z-direction

# Some timestepping parameters
max_Δt = 3600 # maximum allowable timestep 
duration = 100days # The non-dimensional duration of the simulation

# Set the amplitude of the random perturbation (kick)
kick = 0.01  # random velocty perturbation

f = 1e-4 # Coriolis parameter
g = 9.81 # gravitational acceleration
ρ₀ = 1020 # reference density

M2 = 1e-7 # Initial horizontal buoyancy gradient
N2 = 1e-4 # Initial vertical buoyancy gradient

# Set the domain size in non-dimensional coordinates
Lz = 100   # domain size in the vertical (z) direction 
Us=M2*Lz/f  # Set the thermal wind velocity scale
Ri=N2*f^2/M2^2  # Set the Richardson number
Ls = 2*pi*Us/f*sqrt((1+Ri)/(5/2))  # The most unstable mode of BCI from Stone 1970
Lx = 4 * Ls  # size in the x-direction
Ly = 4 * Ls  # size in the x-direction

# Define the background buoyancy
parameters = (N2=N2, M2=M2, Lz=Lz, f=f) # N2 and M2 act as parameters
B(x, y, z, t, p) = p.M2 * y 
U(x, y, z, t, p) = -p.M2 / p.f * (z + p.Lz)
B_field = BackgroundField(B, parameters=parameters)
U_field = BackgroundField(U, parameters=parameters)

# construct a rectilinear grid using an inbuilt Oceananigans function
# Here, the topology parameter sets the style of boundaries in the x, y, and z directions
# 'Bounded' corresponds to wall-bounded directions and 'Flat' corresponds to the dimension that is not considered (here, that is the y direction)
grid = RectilinearGrid(size = (Nx, Ny, Nz), x = (0, Lx), y = (0, Ly), z = (-Lz, 0), topology = (Periodic, Periodic, Bounded))

b_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition(N2),
                                bottom = GradientBoundaryCondition(N2))

# Now, define a 'model' where we specify the grid, advection scheme, bcs, and other settings
model = NonhydrostaticModel(; grid,
              advection = UpwindBiasedFifthOrder(),  # Specify the advection scheme.  Another good choice is WENO() which is more accurate but slower
            timestepper = :RungeKutta3, # Set the timestepping scheme, here 3rd order Runge-Kutta
                tracers = :b,  # Set the name(s) of any tracers, here b is buoyancy and c is a passive tracer (e.g. dye)
               buoyancy = Buoyancy(model=BuoyancyTracer()), # this tells the model that b will act as the buoyancy (and influence momentum) 
               boundary_conditions = (b = b_bcs,),
               background_fields = (b = B_field, u = U_field),
               coriolis = coriolis = FPlane(f = f)
)

# Set initial conditions
# Here, we start with a tanh function for buoyancy and add a random perturbation to the velocity. 
uᵢ(x, y, z) = kick * randn()
vᵢ(x, y, z) = kick * randn() 
wᵢ(x, y, z) = kick * randn()
bᵢ(x, y, z) = N2 * z

# Send the initial conditions to the model to initialize the variables
set!(model, u = uᵢ, v = vᵢ, w = wᵢ, b = bᵢ)

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

ζ = ∂x(v) - ∂y(u) # The vertical vorticity

# Set the name of the output file
filename = "Project7/raw-output/BI_xz"
simulation.output_writers[:xz_slices] =
    JLD2OutputWriter(model, (; u, v, w, b, ζ),
                          filename = filename * ".jld2",
                          indices = (:, 1, :),
                         schedule = TimeInterval(1day),
                            overwrite_existing = true)

filename = "Project7/raw-output/BI_xy"
simulation.output_writers[:xy_slices] =
    JLD2OutputWriter(model, (; u, v, w, b, ζ),
                          filename = filename * ".jld2",
                          indices = (:, :, Nz),
                        schedule = TimeInterval(1day),
                            overwrite_existing = true)

filename = "Project7/raw-output/BI_yz"
simulation.output_writers[:yz_slices] =
    JLD2OutputWriter(model, (; u, v, w, b, ζ),
                         filename = filename * ".jld2",
                           indices = (1, :, :),
                        schedule = TimeInterval(1day),
                            overwrite_existing = true)

nothing # hide

# Now, run the simulation
run!(simulation)

# After the simulation is different, plot the results and save a movie
include("plot_BI_john.jl")
