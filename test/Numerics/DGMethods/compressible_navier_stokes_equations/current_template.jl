include(pwd() * "/test/Numerics/DGMethods/compressible_navier_stokes_equations/boiler_plate.jl")

Ω = Periodic(-2π, 2π)^2 × Interval(-2π, 2π)
grid = DiscretizedDomain(
    Ω,
    elements = (vertical = 1, horizontal = 1),
    polynomialorder = (vertical = 1, horizontal = 1)
)

ρu_bcs = (bottom = Impenetrable(NoSlip()), top = Impenetrable(FreeSlip()))
ρθ_bcs = (bottom = Insulating(), top = TemperatureFlux((state, aux, t) -> 1e-6))
boundary_conditions = (ρθ = ρθ_bcs, ρu = ρu_bcs)

parameters = (
    ϵ  = 0.1, # perturbation size for initial condition
    ρₒ = 1.0, # reference density
    cₛ = 1.0, # linearized sound speed in horizontal
    cᶻ = 1.0, # linearized sound speed in vertical
    α  = 0.0, # thermal expansion coefficient
    g  = 0.0, # gravity
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

ρu₀(x, y, z, p) = sin(x) * sin(y) 
ρv₀(x, y, z, p) = p.ϵ * sin(x) * sin(y) 
ρw₀(x, y, z, p) = 0.0
momentum₀(x,y,z,p) = @SVector [ρu₀(x,y,z,p), ρv₀(x,y,z,p), ρw₀(x,y,z,p)]
initial_conditions = (
    ρ = 1,
    ρu = momentum₀,
    ρθ = (x,y,z,p) -> (0.01 * z),
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