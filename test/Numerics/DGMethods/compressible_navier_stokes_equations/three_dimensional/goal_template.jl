# all the stuff that always needs to be there, using ClimateMachine, ClimateMachine init and so forth
include("boiler_plate.jl") 

# Domain, create tensor product domain [-2π,2π) × [-π,π) × [0,1]
# can also have constructors for Ω = Atmosphere() or Ω = Ocean() etc.
Ω = Periodic(-2π,2π) × Periodic(-π,π) × Interval(0,1) 
# Grid with smart Defaults (like topology selection, float type and so forth)
grid = DiscontinuousSpectralElementGrid(
    domain = Ω,
    elements = (vertical = 4, horizontal = 8)
    polynomialorder = (vertical = 1, horizontal = 5)
)

# Boundary Conditions
# velocity boundary conditions, default free-slip
ρu_bcs = (bottom = NoSlip(), top = FreeSlip())
# temperature, default no-flux
ρθ_bcs = (bottom = Neumann(0.01), top = Flux(1e-6))
# combine
boundary_conditions = (ρθ = ρθ_bcs, ρu = ρu_bcs)

# parameters nobs,
parameters = (
    ϵ = 0.1,
    ρ₀ = 1.0,
    cₛ = 1.0, 
)

# Numerics Nobs
flux = RoeNumericalFlux()

# Dissipation models
ν = IsotropicLaplacian(1e-2)
κ = IsotropicLaplacian(1e-4)
dissipation = (ρu = ν, ρθ = κ)

# Construct the spatial model (implied by balance law)
model = ThreeDimensionalCompressibleNavierStokes(
    grid = grid,
    boundary_conditions = boundary_conditions,
    parameters = parameters,
    dissipation = dissipation,
    flux = flux,
)

# Construct the time evolution
# Initial Conditions, with functions, anonymous functions, scalars, etc.
# can have dependence on parameter struct, p
u₀(x,y,z) = sin(x)*sin(y)
v₀(x, y, z, p) = p.ϵ * sin(x)*sin(y) 
initial_conditions = (
    ρ = 1,
    ρu = (u₀, v₀, 0),
    ρθ = (x,y,z) -> 0.01 * z,
)

# cfl and timestepper
cfl_dt = calculate_cfls(grid, wavespeed = parameters.cₛ)
timestepper = SSPRK22Heuns(dt = cfl_dt)

# Callback vector
callbacks = (
    StateCheck(minute = 1),
    Checkpointer(iteration = 1000),
    MakieVisuals(minute = 10, directory = "/wherever/"),
)

simulation = Simulation(
    model = model,
    inital_conditions = initial_conditions,
    timestepper = timestepper,
    callbacks = callbacks,
)

# and run it
evolve!(simulation)