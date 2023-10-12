# A script to explore Rossby wave propagation and nonlinear behavior using Oceananigans

# Load some packages that we will need
using Oceananigans
using Printf
using Oceananigans.Units
using JLD2

# Set the domain size in dimensional coordinates
Lx = 100kilometers  
Ly = 300kilometers
Lz = 50meters

# Set the grid size
Nx = 32   # number of grid points in x-direction
Ny = 128   # number of grid points in y-direction
Nz = 32

# Construct a rectilinear grid that is periodic in x-direction and bounded in y-direction
grid = RectilinearGrid(size = (Nx, Ny, Nz),
                       x = (0, Lx), y = (0, Ly), z= (-Lz,0),
                       topology = (Periodic, Bounded, Bounded)
)

# Set up a model for Rossby waves
model = NonhydrostaticModel(; grid,
              advection = UpwindBiasedFifthOrder(),   # Specify the advection scheme.  Another good choice is WENO() which is more accurate but slower
            timestepper = :RungeKutta3,   # Set the timestepping scheme, here 3rd order Runge-Kutta
                tracers = :b,
                buoyancy = Buoyancy(model=BuoyancyTracer()),
                coriolis = BetaPlane(rotation_rate = 7.292115e-5, latitude = 45, radius = 6371e3)   # set Coriolis parameter using the Beta-plane approximation 
)

# Set wavenumbers associated with the initial condition
M2 = 1e-6
N2 = 1e-4
f = 7.292115e-5 * sin(pi/4)
kick = 0.05

# Define functions for the initial conditions
uᵢ(x, y, z) = -M2/f*(z+Lz) + kick*randn()
vᵢ(x, y, z) = kick * randn()
wᵢ(x, y, z) = kick * randn()
bᵢ(x, y, z) = N2*z+M2*y # Here, we set the function for c so that it is proportional to the streamfunction associated with (u,v)

# Send the initial conditions to the model to initialize the variables
set!(model, u = uᵢ, v = vᵢ, w = wᵢ, b = bᵢ)

# Create a 'simulation' to run the model for a specified length of time
simulation = Simulation(model, Δt = 10minutes, stop_iteration = 200)

# Add callback that prints progress message during simulation
progress(sim) = @info string("Iter: ", iteration(sim),
                             ", time: ", prettytime(sim))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

# Save output from the simulation
filename = "rossbywave"

u, v, w = model.velocities
b = model.tracers.b
ζ = ∂x(v)-∂y(u)

simulation.output_writers[:jld2] = JLD2OutputWriter(model, (; u, v, w, b, ζ),
                                                    schedule = IterationInterval(10),
                                                    filename = filename * ".jld2",
                                                    indices = (:,:,1),
                                                    overwrite_existing = true
)

# Run the simulation                                                  
run!(simulation)

# Make a plot of u at y=Ly/2 and save a movie
include("plot_rossbywave_edit.jl")
