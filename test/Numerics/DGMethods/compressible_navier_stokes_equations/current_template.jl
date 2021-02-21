include(pwd() * "/test/Numerics/DGMethods/compressible_navier_stokes_equations/boiler_plate.jl")

Ω = Periodic(-2π,2π) × Periodic(-2π,2π) × Periodic(-2π,2π)
grid = DiscretizedDomain(
    Ω,
    elements = (vertical = 1, horizontal = 1),
    polynomialorder = (vertical = 1, horizontal = 1)
)

ρu_bcs = (bottom = Impenetrable(NoSlip()), top = Impenetrable(FreeSlip()))
ρθ_bcs = (bottom = TemperatureFlux(0.00), top = TemperatureFlux(1e-6))
boundary_conditions = (ρθ = ρθ_bcs, ρu = ρu_bcs)

parameters = (
    ϵ = 0.1,
    ρₒ = 1.0,
    cₛ = 1.0, 
    cᶻ = 1.0,
    α = 0.0,
    g = 0.0,
)

flux = RoeNumericalFlux()

ν = Laplacian(1e-2)
κ = Laplacian(1e-4)
dissipation = (ρu = ν, ρθ = κ)

model = SpatialModel(
    balancelaw = ThreeDimensionalCompressibleNavierStokesEquations,
    physics = (;dissipation, ),
    numerics = (;flux),
    grid = grid,
    boundaryconditions = boundary_conditions,
    parameters = parameters
)

ρu₀(x, y, z, p) = state.ρ * sin(x) * sin(y) 
ρv₀(x, y, z, p) = state.ρ * p.ϵ * sin(x) * sin(y) 
initial_conditions = (
    ρ = 1,
    ρu = (ρu₀, ρv₀, 0),
    ρθ = (x,y,z,p) -> state.ρ * (0.01 * z),
)

Δt = calculate_dt(grid, wavespeed = parameters.cₛ, cfl = 0.1)
timestepper = TimeStepper(method = SSPRK22Heuns, timestep = Δt)

simulationtime = (0, 10)

callbacks = (
    Default(),
)

simulation = Simulation(
    model = model,
    initialconditions = initial_conditions,
    timestepper = timestepper,
    callbacks = callbacks,
    simulationtime = simulationtime,
)