include("boilerplate.jl")
include("two_dimensional/TwoDimensionalCompressibleNavierStokesEquations.jl")
include("three_dimensional/ThreeDimensionalCompressibleNavierStokesEquations.jl")
ClimateMachine.init(disable_gpu = true)

########
# Setup physical and numerical domains
########
Ω = Periodic(-100, 100)^2 × Interval(-100, 0)
grid = DiscretizedDomain(
    Ω,
    elements = (vertical = 1, horizontal = 1),
    polynomialorder = (vertical = 1, horizontal = 1),
)

########
# Define physical parameters and parameterizations
########
FT = eltype(grid.numerical.vgeo)
Lx, Ly, Lz = length(Ω)
parameters = (
    ρₒ = 1.0,  # reference density
    c = 2,     # Rusanov wavespeed
    g = 10,    # gravity
    cₛ = 1.0,  # horizontal linearized sound speed
    cᶻ = 1.0,  # vertical linearized sound speed
    α = 2e-4,  # should probably multiply by cᶻ^2 due to hydrostatic balance
    ν = 1e-2,  # kg / (m s)
    κ = 1e-4,  # kg / (m s)
    Lx = Lx,   # domain length in horizontal
    Ly = Ly,   # domain length in horizontal
    Lz = Lz,   # domain length in vertical
)

dissipation = ConstantViscosity{Float64}(
    μ = 0, 
    ν = parameters.ν, 
    κ = parameters.κ
)

physics = FluidPhysics(;
    advection = NonLinearAdvectionTerm(),
    dissipation = dissipation,
    coriolis = nothing,
    buoyancy = Buoyancy{FT}(α = parameters.α, g = parameters.g),
)

########
# Define boundary conditions and numerical fluxes
########
ρu_bc = Impenetrable(FreeSlip())
ρθ_bc = Insulating()
top_ρθ_bc = TemperatureFlux((state, aux, t) -> 1e-5)
ρu_bcs = (south = ρu_bc, north = ρu_bc)
ρθ_bcs = (south = ρθ_bc, north = ρθ_bc)
BC = (ρθ = ρθ_bcs, ρu = ρu_bcs)

flux = RoeNumericalFlux()
overintegration = 0

numerics = (; flux, overintegration)

########
# Define initial conditions
########

ρ₀(x, y, z, p) = p.ρₒ * ( 1 +  ( p.α * p.g / p.cᶻ^2) * z^2 / (2 * p.Lz))
ρu₀(x, y, z, p) = ρ₀(x, y, z, p) * p.ϵ * sin(2π*y/p.Ly)*sin(2π*z/p.Lz) 
ρv₀(x, y, z, p) = ρ₀(x, y, z, p) * p.ϵ * sin(2π*x/p.Lx)*sin(2π*z/p.Lz)
ρw₀(x, y, z, p) = ρ₀(x, y, z, p) * 0.0
ρθ₀(x, y, z, p) = ρ₀(x, y, z, p) * (z / p.Lz )

ρu⃗₀(x, y, z, p) =
    @SVector [ρu₀(x, y, z, p), ρv₀(x, y, z, p), ρw₀(x, y, z, p)]
initial_conditions = (ρ = ρ₀, ρu = ρu⃗₀, ρθ = ρθ₀)

########
# Define timestepping parameters
########
days = 86400.0
start_time = 0
end_time = 0.5days
method = SSPRK22Heuns

Δt = calculate_dt(grid, wavespeed = parameters.cᶻ, cfl = 0.2)

########
# Define callbacks
########
iteration = floor(Int, end_time / Δt)
jldcallback = JLD2State(iteration = iteration, filepath = "convection.jld2")
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
println("The amount of time for the simulation was ", toc -tic, " seconds")

##
visualize(simulation)
