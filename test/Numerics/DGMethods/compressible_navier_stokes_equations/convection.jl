include("boilerplate.jl")
include("two_dimensional/TwoDimensionalCompressibleNavierStokesEquations.jl")
include("three_dimensional/ThreeDimensionalCompressibleNavierStokesEquations.jl")
ClimateMachine.init()

########
# Setup physical and numerical domains
########
Ω = Periodic(-100, 100)^2 × Interval(-100, 0)
grid = DiscretizedDomain(
    Ω,
    elements = (vertical = 1, horizontal = 8),
    polynomialorder = (vertical = 1, horizontal = 3),
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
    cₛ = 2,
    cᶻ = 3,
)

physics = FluidPhysics(;
    advection = NonLinearAdvectionTerm(),
    dissipation = ConstantViscosity{Float64}(μ = 0, ν = 0, κ = 0),
    coriolis = nothing,
    buoyancy = Buoyancy{FT}(α = parameters.α, g = parameters.g),
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
end_time = 1.0
method = SSPRK22Heuns

Δt = calculate_dt(grid, wavespeed = sqrt(parameters.g), cfl = 0.3)

########
# Define callbacks
########

jldcallback = JLD2State(iteration = 100, filepath = "test.jld2")
callbacks = (Info(), StateCheck(10), jldcallback)
callbacks = (Info(),)

########
# Create the things
########
model = SpatialModel(
    balance_law = Fluid3D(),
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
tic = Base.time()
evolve!(simulation, model)
toc = Base.time()
println("The amount of time for the simulation was ", toc -tic)

##
visualize(simulation)
