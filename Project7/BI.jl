# Non-dimensionalised Eady model

# Load in Oceannanigans and some standard libraries
using Printf
using Oceananigans
using Oceananigans.Units

# We have non-dimensionalised so that f = 1 is the Coriolis parameter and H = 1 is the depth of the mixed layer
# Set the two parameters for the problem
Ri = 100       # Richardson number, N²f²/M⁴
α = 10         # M²/f²
# Calculate some other non-dimensional parameters
β = Ri * α^2   # N²/f²
Ro = Ri^(-0.5) # (initial) Rossby number

# Print the non-dimensional numbers
@printf("\nThe non-dimensional numbers are:\n")
@printf("Ri = %.2e\n", Ri)
@printf("α = %.2e\n", α)
@printf("β = %.2e\n", β)
@printf("Ro = %.2e\n\n", Ro)

# Now calculate dimensionalised parameters in SI (under the assumption f = 1e-4 and H = 50) to help with interpretation
f_dim = 1e-4
H = 50
M²_dim = α * f_dim^2                        # The horizontal buoyancy gradient
N²_dim = β * f_dim^2                        # The vertical buoyancy gradient
U_dim = α * f_dim * H                       # The initial velocity magnitude
L₀ = β^0.5 * H                              # The Rossby deformation radius
L_dim = 2*pi*U_dim/f_dim * ((1+Ri)*2/5)^0.5 # The wavelength of the most unstable mode (Stone, 1966)
τ_dim = (54/5)^0.5 * (1+Ri)^0.5 / f_dim     # The Timescale of instability growth (Stone, 1966)
# Print the dimensionalised parameters
@printf("In SI units with f = %.2e and H = %i:\n", f_dim, H)
@printf("M² = %.2e\n", M²_dim)
@printf("N² = %.2e\n", N²_dim)
@printf("U = %.2e\n", U_dim)
@printf("L₀ = %.2fkm\n", L₀/1kilometer)
@printf("L = %.2fkm\n", L_dim/1kilometer)
@printf("τ = %.2fh\n", τ_dim/1hour)

# Set the grid size
Nx = 128
Ny = 128
Nz = 16
# Set the domain size
Lx = L_dim/L₀ * 4
Ly = Lx
Lz = 1
# Calculate the non-dimensional growth rate (Stone, 1966)
τ = τ_dim * f_dim
# Use the growth rate to set the simulation timestep and duration
max_Δt = τ / 25
duration = τ * 15

# Define the background buoyancy and velocity fields
B(x,y,z,t) = α * y
U(x,y,z,t) = -α * (z+Lz)
B_field = BackgroundField(B)
U_field = BackgroundField(U)

# Construct a rectilinear grid which is periodic in the x and y directions
grid = RectilinearGrid(size = (Nx, Ny, Nz), x = (0, Lx), y = (0, Ly), z = (-Lz, 0), topology = (Periodic, Periodic, Bounded))

# Set the only non-default boundary conditions, so that the initial conditions satisfy the boundary conditions
b_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition(β),
                                bottom = GradientBoundaryCondition(β))

# Now, define a 'model' where we specify the grid, advection scheme, bcs, and other settings
model = NonhydrostaticModel(; grid,
                advection = UpwindBiasedFifthOrder(),  # Specify the advection scheme.  Another good choice is WENO() which is more accurate but slower
                timestepper = :RungeKutta3, # Set the timestepping scheme, here 3rd order Runge-Kutta
                tracers = :b,  # Set the name(s) of any tracers, here b is buoyancy and c is a passive tracer (e.g. dye)
                buoyancy = Buoyancy(model=BuoyancyTracer()), # This tells the model that b will act as the buoyancy (and influence momentum) 
                boundary_conditions = (b = b_bcs,), # Not sure how the ,) works, but it does
                background_fields = (b = B_field, u = U_field),
                coriolis = coriolis = FPlane(f = 1)
)

uᵢ(x, y, z) = α * 0.05 * randn()
vᵢ(x, y, z) = α * 0.05 * randn() 
wᵢ(x, y, z) = α * 0.05 * randn()
bᵢ(x, y, z) = β * z

# Now, we create a 'simulation' to run the model for a specified length of time
simulation = Simulation(model, Δt = max_Δt/5, stop_time = duration)

# ### The `TimeStepWizard`
#
# The TimeStepWizard manages the time-step adaptively, keeping the
# Courant-Freidrichs-Lewy (CFL) number close to `1.0` while ensuring
# the time-step does not increase beyond the maximum allowable value
wizard = TimeStepWizard(cfl = 0.5, max_change = 1.1, max_Δt = max_Δt)
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
# Field has to be used for outputs which get processed because otherwise a bug occurs - Jago moment

u = Field(model.velocities.u + model.background_fields.velocities.u) # unpack velocity `Field`s
v = model.velocities.v
w = model.velocities.w
b = Field(model.tracers.b + model.background_fields.tracers.b) # extract the buoyancy and add the background field
ζ = Field(∂x(v) - ∂y(u))    # The vertical vorticity
δ = Field(∂x(u)+∂y(v))      # The horizontal divergence

# Set the name of the JLD2 output file
filename = "Project7/raw-output/BI_xz"
# Output the slice y = 0
simulation.output_writers[:xz_slices] =
    JLD2OutputWriter(model, (; u, v, w, b, ζ, δ),
                            filename = filename * ".jld2",
                            indices = (:, 1, :),
                            schedule = TimeInterval(1day),
                            overwrite_existing = true)

# Set the name of the JLD2 output file
filename = "Project7/raw-output/BI_xy"
# Output the slice z = 0 (top surface)
simulation.output_writers[:xy_slices] =
    JLD2OutputWriter(model, (; u, v, w, b, ζ, δ),
                            filename = filename * ".jld2",
                            indices = (:, :, Nz),
                            schedule = TimeInterval(1day),
                            overwrite_existing = true)

# Set the name of the JLD2 output file
filename = "Project7/raw-output/BI_yz"
# Output the slice x = 0
simulation.output_writers[:yz_slices] =
    JLD2OutputWriter(model, (; u, v, w, b, ζ, δ),
                            filename = filename * ".jld2",
                            indices = (1, :, :),
                            schedule = TimeInterval(τ/20),
                            overwrite_existing = true)

nothing # Hide

# Now, run the simulation
run!(simulation)

# After the simulation is done, plot the results and save a movie
include("plot_BI_john.jl")

@printf("\nEnd")