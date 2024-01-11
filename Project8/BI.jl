# Eady model

# Load some standard libraries that we will need
using Printf
using Oceananigans
using Oceananigans.Units
using Oceananigans.AbstractOperations

include("plot_BI.jl")

function run_BI_sim(Ri, α, label)

# Ri = N²f²/M⁴ is the linear instability parameter
# α = M²/N² is the slope of the isopycnals
# label is a suffix for saving data from the simulations

# Set the grid size
Nx = 128  # number of gridpoints in the x-direction
Ny = 128  # number of gridpoints in the x-direction
Nz = 16   # number of gridpoints in the z-direction

# Set the two dimensional parameters
H = 50    # Depth of mixed layer
f = 1e-4  # Coriolis parameter

# Calculate the other non-dimensional parameters
s = α^2 / Ri        # N²/f²
λ = α / Ri          # M²/f²
Ro₀ = Ri ^ -0.5     # The initial Rossby number

# Calculate the other dimensional parameters and scales
M² = λ * f^2                # Horizontal buoyancy gradient
N² = s * f^2                # Vertical buoyancy gradient
L = λ * H * (1+Ri) ^ 0.5    # Horizontal lengthscale (Rossby deformation radius)
U = λ * f * H               # Velocity scale
T = (1+Ri) ^ 0.5 / f        # Timescale of growth (and expected timescale of later flow, if Ro ∼ Ri^(-0.5)

@info M², N²

# Set the domain size
Lx = 4 * 2*pi * L * 0.4^0.5 # Zonal extent, set to 4 wavelengths of the most unstable mode
Ly = (1+5^0.5)/2 * Lx       # Meridional extent, chosen to be an irrational multiple of Lx to avoid resonance
Lz = H                      # Vertical extent

# Now come up with a useful viscosity to use
#=Lx_scale = 4 * 2*pi * s^0.5 * H * 0.4^0.5   # A horizontal lengthscale which only depends on the non-dimensional stratification (not Ri)
ν_h = f * Lx_scale^2 / Nx^2                 # Horizontal viscosity
ν_v = f * H^2 / Nz^2                        # Vertical viscosity=#
ν_h = 1e-4
ν_v = 1e-6
diff_h = HorizontalScalarDiffusivity(ν = ν_h, κ = ν_h)
diff_v = HorizontalScalarDiffusivity(ν = ν_v, κ = ν_v)
@info ν_h
@info ν_v

# Set some timestepping parameters
max_Δt = minimum([T/10, 0.5 * pi / (N²^0.5)])
duration = T * 40

@info U, L, M², N²

@printf("Simulation will last %s\n", prettytime(duration))

# Set relative amplitude for random velocity perturbation
kick = 0.02 * H * f

# Define the background fields
B₀(x, y, z, t) = M² * y + N² * z    # Buoyancy
U₀(x, y, z, t) = -M²/f * (z + Lz)   # Zonal velocity
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
              buoyancy = Buoyancy(model = BuoyancyTracer()), # this tells the model that b will act as the buoyancy (and influence momentum)
              background_fields = (b = B_field, u = U_field),
              coriolis = coriolis = FPlane(f = f),
              closure = (diff_h, diff_v)
)

# Set initial conditions, a random velocity perturbation
uᵢ(x, y, z) = kick * randn()
vᵢ(x, y, z) = kick * randn() 
wᵢ(x, y, z) = kick * randn()
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
wizard = TimeStepWizard(cfl = 0.5, max_change = 1.1, max_Δt = max_Δt)

# Still some numerical noise at CFL 0.1 for Ri = 10⁴, but none for CFL = 0.05

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

u = Field(model.velocities.u + model.background_fields.velocities.u)    # Unpack velocity `Field`s
v = Field(model.velocities.v)
w = Field(model.velocities.w)
b = Field(model.tracers.b + model.background_fields.tracers.b)          # Extract the buoyancy and add the background field
u_pert = Field(model.velocities.u)
b_pert = Field(model.tracers.b)


ζ = Field(∂x(v) - ∂y(u))    # The vertical vorticity
δ = Field(∂x(u) + ∂y(v))    # The horizontal divergence
ℬ = Field(w * b_pert)       # The buoyancy flux

# Output the slice y = 0
filename = "Project8/raw-output/BI_xz" * label
simulation.output_writers[:xz_slices] =
    JLD2OutputWriter(model, (; u, v, w, b, ζ, δ),
                            filename = filename * ".jld2",
                            indices = (:, 1, :),
                            schedule = TimeInterval(T/20),
                            overwrite_existing = true)

# Output the slice z = 0
filename = "Project8/raw-output/BI_xy" * label
simulation.output_writers[:xy_slices] =
    JLD2OutputWriter(model, (; u, v, w, b, ζ, δ),
                            filename = filename * ".jld2",
                            indices = (:, :, Nz),
                            schedule = TimeInterval(T/20),
                            overwrite_existing = true)

# Output the slice x = 0
filename = "Project8/raw-output/BI_yz" * label
simulation.output_writers[:yz_slices] =
    JLD2OutputWriter(model, (; u, v, w, b, ζ, δ),
                            filename = filename * ".jld2",
                            indices = (1, :, :),
                            schedule = TimeInterval(T/20),
                            overwrite_existing = true)

# Output the total vertical buoyancy flux
filename = "Project8/raw-output/BI_ℬ" * label
simulation.output_writers[:ℬ_flux] =
    JLD2OutputWriter(model, (; ℬ),
                            filename = filename * ".jld2",
                            indices = (:, :, Int64(round((Nz+1)/2))),
                            schedule = TimeInterval(T/20),
                            overwrite_existing = true)


nothing # hide

# Now, run the simulation
run!(simulation)

# After the simulation is done, plot the results and save a movie

BI_plot(label)

end