#!/usr/bin/env julia --project
include("../boilerplate.jl")
include("TwoDimensionalCompressibleNavierStokesEquations.jl")

ClimateMachine.init()

#################
# RUN THE TESTS #
#################

########
# Setup physical and numerical domains
########
Ω = Periodic(-2π, 2π) × Interval(-2π, 2π)
grid = DiscretizedDomain(Ω, elements = 16, polynomialorder = 3)

########
# Define physical parameters and parameterizations
########
FT = eltype(grid.numerical.vgeo)

parameters = (
    ϵ  = 0.1,  # perturbation size for initial condition
    l  = 0.5,  # Gaussian width
    k  = 0.5,  # Sinusoidal wavenumber
    ρₒ = 1.0, # reference density
    c  = 2,
    g  = 10,
)

dissipation = ConstantViscosity{FT}(μ = 0, ν = 0, κ = 0)

physics = FluidPhysics(;
    advection = NonLinearAdvectionTerm(),
    dissipation = dissipation,
    coriolis = nothing,
    buoyancy = nothing,
)

########
# Define boundary conditions and numerical fluxes
########
ρu_bc = Impenetrable(NoSlip())
ρθ_bc = Insulating()
ρu_bcs = (south = ρu_bc, north = ρu_bc)
ρθ_bcs = (south = ρθ_bc, north = ρθ_bc)
BC = (ρθ = ρθ_bcs, ρu = ρu_bcs)

flux = RoeNumericalFlux()
overintegration = 1

numerics = (; flux, overintegration)

########
# Define initial conditions
########

# The Bickley jet
U₀(x, y, z, p) = cosh(y)^(-2)

# Slightly off-center vortical perturbations
Ψ₀(x, y, z, p) =
    exp(-(y + p.l / 10)^2 / (2 * (p.l^2))) * cos(p.k * x) * cos(p.k * y)

# Vortical velocity fields (ũ, ṽ) = (-∂ʸ, +∂ˣ) ψ̃
u₀(x, y, z, p) = Ψ₀(x, y, z, p) * (p.k * tan(p.k * y) + y / (p.l^2))
v₀(x, y, z, p) = -Ψ₀(x, y, z, p) * p.k * tan(p.k * x)

ρ₀(x, y, z, p) = p.ρₒ
ρu₀(x, y, z, p) = ρ₀(x, y, z, p) * (U₀(x, y, z, p) + p.ϵ * u₀(x, y, z, p))
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
end_time = 80.0
Δt = 0.02
method = SSPRK22Heuns

########
# Define callbacks
########

callbacks = () # (Info(), )

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
    simulation_time = (start_time, end_time),
)

########
# Run the model
########
evolve!(simulation, model)

########
# Visualize the Model
########

scene = visualize(simulation, statenames = ["ρ", "ρu", "ρv", "ρθ"])

#######
# Record Interaction
#######
seconds = 10
fps = 30
frames = round(Int, fps * seconds )
record(scene, pwd() * "/example.mp4"; framerate = fps) do io
    for i = 1:frames
        sleep(1/fps)
        recordframe!(io)
    end
end

