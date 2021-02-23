include("boilerplate.jl")

########
# Setup physical and numerical domains
########
Ω = Periodic(-2π, 2π) × Periodic(-2π, 2π) × Interval(-2π, 2π)
grid = DiscretizedDomain(
    Ω,
    elements = (vertical = 2, horizontal = 2),
    polynomialorder = (vertical = 1, horizontal = 1),
)

########
# Define boundary conditions and numerical fluxes
########
ρu_bcs = (bottom = Impenetrable(NoSlip()), top = Impenetrable(FreeSlip()))
ρθ_bcs = (bottom = Insulating(), top = TemperatureFlux((state, aux, t) -> 1e-6))
boundary_conditions = (ρθ = ρθ_bcs, ρu = ρu_bcs)

flux = RoeNumericalFlux()

numerics = (; flux)

########
# Define physical parameters and parameterizations
########
FT = eltype(grid.numerical.vgeo)

parameters = (
    ϵ = 0.1, # perturbation size for initial condition
    ρₒ = 1.0, # reference density
    cₛ = 1.0, # linearized sound speed in horizontal
    cᶻ = 1.0, # linearized sound speed in vertical 
)

buoyancy = Buoyancy{FT}(α = 0.0, g = 0.0)
dissipation = ConstantVelocity{FT}(μ = 0, ν = 1e-6, κ = 1e-12)

physics = FluidPhysics(;
    advection = NonLinearAdvectionTerm(),
    dissipation = dissipation,
    coriolis = nothing,
    buoyancy = nothing,
)

########
# Define initial conditions
########
ρ₀(x, y, z, p) = p.ρₒ
ρu₀(x, y, z, p) = ρ(x, y, z, p) * sin(x) * sin(y)
ρv₀(x, y, z, p) = ρ(x, y, z, p) * p.ϵ * sin(x) * sin(y)
ρw₀(x, y, z, p) = ρ(x, y, z, p) * 0.0
ρθ₀(x, y, z, p) = ρ(x, y, z, p) * (0.01 * z)

ρu₀(x, y, z, p) = @SVector [ρu₀(x, y, z, p), ρv₀(x, y, z, p), ρw₀(x, y, z, p)]
initial_conditions = (ρ = ρ₀, ρu = ρu₀, ρθ = ρθ₀)

########
# Define timestepping parameters
########
start_time = 0
end_time = 0
method = SSPRK22Heuns

Δt = calculate_dt(grid, wavespeed = parameters.cₛ, cfl = 0.1)

########
# Define callbacks
########

callbacks = (Info())

########
# Create the things
########
model = SpatialModel(
    balance_law = CNSE3D,
    physics = physics,
    numerics = numerics,
    grid = grid,
    boundary_conditions = boundary_conditions,
    parameters = parameters,
)

timestepper = TimeStepper(method = method, timestep = Δt)

simulation = Simulation(
    model = model,
    initial_conditions = initial_conditions,
    timestepper = timestepper,
    callbacks = callbacks,
    simulation_time = (start_time, end_time),
)

########
# Run the model
########

evolve!(simulation)
