include("boilerplate.jl")
include("two_dimensional/TwoDimensionalCompressibleNavierStokesEquations.jl")

ClimateMachine.init()

########
# Setup physical and numerical domains
########
Ω = Periodic(-2π, 2π) × Periodic(-2π, 2π)
grid = DiscretizedDomain(
    Ω,
    elements = (vertical = 8, horizontal = 8),
    polynomialorder = (vertical = 3, horizontal = 3),
)

########
# Define physical parameters and parameterizations
########
FT = eltype(grid.numerical.vgeo)

parameters = (
    ϵ = 0.1,  # perturbation size for initial condition
    l = 0.5, # Gaussian width
    k = 0.5, # Sinusoidal wavenumber
    ρₒ = 1.0, # reference density
    c = 2,
    g = 10,
)

physics = FluidPhysics(;
    advection = NonLinearAdvectionTerm(),
    dissipation = ConstantViscosity{Float64}(μ = 0, ν = 0, κ = 0),
    coriolis = nothing,
    buoyancy = nothing,
)


########
# Define boundary conditions and numerical fluxes
########
ρu_bc = Impenetrable(FreeSlip())
ρθ_bc = Insulating()
ρu_bcs = (south = ρu_bc, north = ρu_bc)
ρθ_bcs = (south = ρθ_bc, north = ρθ_bc)
BC = (ρθ = ρθ_bcs, ρu = ρu_bcs)

flux = RoeNumericalFlux()

numerics = (; flux)

########
# Define initial conditions
########
U₀(x, y, z, p) = cosh(y)^(-2)

# Slightly off-center vortical perturbations
Ψ₀(x, y, z, p) =
    exp(-(y + p.l / 10)^2 / (2 * (p.l^2))) * cos(p.k * x) * cos(p.k * y)

# Vortical velocity fields (ũ, ṽ) = (-∂ʸ, +∂ˣ) ψ̃
u₀(x, y, z, p) = Ψ₀(x, y, z, p) * (p.k * tan(p.k * y) + y / (p.l^2))
v₀(x, y, z, p) = -Ψ₀(x, y, z, p) * p.k * tan(p.k * x)

ρ₀(x, y, z, p) = p.ρₒ
ρu₀(x, y, z, p) = ρ₀(x, y, z, p) * (p.ϵ * u₀(x, y, z, p) + U₀(x, y, z, p))
ρv₀(x, y, z, p) = ρ₀(x, y, z, p) * p.ϵ * v₀(x, y, z, p)
ρw₀(x, y, z, p) = ρ₀(x, y, z, p) * 0.0
ρθ₀(x, y, z, p) = ρ₀(x, y, z, p) * sin(p.k * y)

ρu⃗₀(x, y, z, p) =
    @SVector [ρu₀(x, y, z, p), ρv₀(x, y, z, p), ρw₀(x, y, z, p)]
initial_conditions = (ρ = ρ₀, ρu = ρu⃗₀, ρθ = ρθ₀)

########
# Define timestepping parameters
########
start_time = 0
end_time = 50.0
method = SSPRK22Heuns

Δt = calculate_dt(grid, wavespeed = sqrt(parameters.g), cfl = 0.3)

########
# Define callbacks
########

callbacks = (Info(), StateCheck(10))

########
# Create the things
########
model = SpatialModel(
    balance_law = Fluid2D(),
    physics = physics,
    numerics = numerics,
    grid = grid,
    boundary_conditions = BC,
    parameters = parameters,
)

timestepper = TimeStepper(method = method, timestep = Δt)

simulation = Simulation(
    model = model,
    initial_conditions = initial_conditions,
    timestepper = timestepper,
    callbacks = callbacks,
    time = (; start = start_time, finish = end_time),
)

########
# Run the model
########

evolve!(simulation, model)
