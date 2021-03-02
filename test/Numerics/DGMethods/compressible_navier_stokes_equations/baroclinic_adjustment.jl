include("boilerplate.jl")
include("two_dimensional/TwoDimensionalCompressibleNavierStokesEquations.jl")
include("three_dimensional/ThreeDimensionalCompressibleNavierStokesEquations.jl")
ClimateMachine.init(disable_gpu = false)
# include("test/Numerics/DGMethods/compressible_navier_stokes_equations/convection.jl")

########
# Setup physical and numerical domains
########
Ω = Periodic(0, 250e3) × Interval(0, 500e3) × Interval(-1e3, 0)

overintegration = 1
grid = DiscretizedDomain(
    Ω,
    elements = (32, 64, 4),
    polynomialorder = (vertical = 1+overintegration, horizontal = 3+overintegration),
)

########
# Define physical parameters and parameterizations
########
FT = eltype(grid.numerical.vgeo)
Lx, Ly, Lz = length(Ω)
scale = 1.0 # Lz / Lx
parameters = (
    ρₒ = 1.0,  # reference density
    ϵ = 0.0,  # perturbation velocity amplitude
    c = 2,     # Rusanov wavespeed
    g = 1.0,    # gravity
    cₛ = 10.0,  # horizontal linearized sound speed
    cᶻ = 1.0,  # vertical linearized sound speed
    α = 2e-4,  # should probably multiply by cᶻ^2 due to hydrostatic balance
    ν = 2e1,  # kg / (m s)
    κ = 2e1,  # kg / (m s)
    Lx = Lx,   # domain length in horizontal
    Ly = Ly,   # domain length in horizontal
    Lz = Lz,   # domain length in vertical
)

dissipation = ConstantViscosity{Float64}(
    μ = 0, 
    ν = parameters.ν, 
    κ = parameters.κ
)

coriolis = fPlaneCoriolis{Float64}(
        fₒ = 1e-4, # Hz
        β = 0.0, # Hz/m
) 

physics = FluidPhysics(;
    advection = NonLinearAdvectionTerm(),
    dissipation = dissipation,
    coriolis = coriolis,
    buoyancy = Buoyancy{FT}(α = parameters.α, g = parameters.g),
)

########
# Define boundary conditions and numerical fluxes
########
ρu_bc = Impenetrable(FreeSlip())
ρθ_bc = Insulating()
top_ρθ_bc = Insulating() #+takes away heat from top
ρu_bcs = (bottom = ρu_bc, top = ρu_bc)
ρθ_bcs = (bottom = ρθ_bc, top = top_ρθ_bc)
BC = (ρθ = ρθ_bcs, ρu = ρu_bcs)

flux = RoeNumericalFlux()

numerics = (; flux, overintegration)

########
# Define initial conditions
########

ρ₀(x, y, z, p) = p.ρₒ * (1 + p.α * p.g * 2 * -z/p.Lz) #( 1 +  ( p.α * p.g / p.cᶻ^2) * (z^2 / p.Lz + (z + p.Lz) * (4e-5*min(max(0, y-225e3), 50e3) - 2 ) ))
ρu₀(x, y, z, p) = ρ₀(x, y, z, p) * p.ϵ * sin(2π*y/p.Ly)*sin(2π*z/p.Lz) 
ρv₀(x, y, z, p) = ρ₀(x, y, z, p) * p.ϵ * sin(2π*x/p.Lx)*sin(2π*z/p.Lz)
ρw₀(x, y, z, p) = ρ₀(x, y, z, p) * p.ϵ * sin(8π*x/p.Ly)*sin(8π*y/p.Ly)
ρθ₀(x, y, z, p) = ρ₀(x, y, z, p) * (-2.0 + 4e-5*min(max(0, y-225e3), 50e3) + 2*z/p.Lz + 0.0001*rand())

ρu⃗₀(x, y, z, p) =
    @SVector [ρu₀(x, y, z, p), ρv₀(x, y, z, p), ρw₀(x, y, z, p)]
initial_conditions = (ρ = ρ₀, ρu = ρu⃗₀, ρθ = ρθ₀)

########
# Define timestepping parameters
########
days = 86400.0
start_time = 0
end_time = 30days
method = SSPRK22Heuns

Δt = calculate_dt(grid, wavespeed = parameters.cᶻ, viscocity = parameters.ν, cfl = 0.2)

########
# Define callbacks
########
iteration = floor(Int, end_time / Δt / 120)
jldcallback = JLD2State(iteration = iteration, filepath = "baroclinic_states.jld2")
callbacks = (Info(), StateCheck(10), jldcallback)
# callbacks = (Info(), StateCheck(10))

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
ρθ = Array(simulation.state.ρθ)
cpu_grid = DiscretizedDomain(
    grid.domain;
    elements = grid.resolution.elements,
    polynomialorder = grid.resolution.polynomialorder,
    array = Array,
)

@save "baroclinic_adjustment.jld2" ρθ cpu_grid
