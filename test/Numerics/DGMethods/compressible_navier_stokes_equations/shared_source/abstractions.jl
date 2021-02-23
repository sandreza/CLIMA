#######
# useful concepts for dispatch
#######

"""
Advection terms

right now really only non-linear or ::Nothing
"""
abstract type AdvectionTerm end
struct NonLinearAdvectionTerm <: AdvectionTerm end

"""
Turbulence Closures

ways to handle drag and diffusion and such
"""
abstract type TurbulenceClosure end

struct LinearDrag{T} <: TurbulenceClosure
    λ::T
end

struct ConstantViscosity{T} <: TurbulenceClosure
    μ::T
    ν::T
    κ::T
    function ConstantViscosity{T}(;
        μ = T(1e-6),   # m²/s
        ν = T(1e-6),   # m²/s
        κ = T(1e-6),   # m²/s
    ) where {T <: AbstractFloat}
        return new{T}(μ, ν, κ)
    end
end

"""
Forcings

ways to add body terms and sources
"""
abstract type Forcing end
abstract type CoriolisForce <: Forcing end

struct fPlaneCoriolis{T} <: CoriolisForce
    fₒ::T
    β::T
    function fPlaneCoriolis{T}(;
        fₒ = T(1e-4), # Hz
        β = T(1e-11), # Hz/m
    ) where {T <: AbstractFloat}
        return new{T}(fₒ, β)
    end
end

struct WindStress{T} <: Forcing
    τₒ::T
    function WindStress{T}(; τₒ = T(1e-4)) where {T <: AbstractFloat}
        return new{T}(τₒ)
    end
end

struct Buoyancy{T} <: Forcing
    α::T # 1/K
    g::T # m/s²
    function Buoyancy{T}(; α = T(2e-4), g = T(10)) where {T <: AbstractFloat}
        return new{T}(α, g)
    end
end

"""
Grouping structs
"""
abstract type AbstractModel end

Base.@kwdef struct SpatialModel{𝒜, ℬ, 𝒞, 𝒟, ℰ, ℱ} <: AbstractModel
    balance_law::𝒜
    physics::ℬ
    numerics::𝒞
    grid::𝒟
    boundary_conditions::ℰ
    parameters::ℱ
end

abstract type ModelPhysics end

Base.@kwdef struct FluidPhysics{A, D, C, B} <: ModelPhysics
    advection::A = NonLinearAdvectionTerm()
    dissipation::D = nothing
    coriolis::C = nothing
    buoyancy::B = nothing
end

abstract type AbstractSimulation end

Base.@kwdef struct Simulation{𝒜, ℬ, 𝒞, 𝒟, ℰ} <: AbstractSimulation
    model::𝒜
    state::ℬ
    timestepper::𝒞
    initial_conditions::𝒟
    callbacks::ℰ
    simulation_time::ℱ
end

function Simulation(;
    model = nothing,
    timestepper = nothing,
    initial_conditions = nothing,
    callbacks = nothing,
    simulation_time = nothing,
)
    model = DGModel(model)

    FT = eltype(simulation.model.grid.numerical.vgeo)
    state = init_ode_state(dg, FT(0); init_on_cpu = true)

    return Simulation(;
        model,
        state,
        timestepper,
        initial_conditions,
        callbacks,
        simulation_time,
    )
end

function DGModel(model::SpatialModel{BL}) where {BL <: CNSE3D}
    params = model.parameters
    physics = model.physics

    Lˣ, Lʸ, Lᶻ = length(model.grid.domain)
    boundary_conditions = get_boundary_conditions(model)

    model = model.balance_law(
        (Lˣ, Lʸ, Lᶻ),
        physics.advection,
        physics.dissipation,
        physics.coriolis,
        physics.buoyancy,
        boundary_conditions,
        ρₒ = params.ρₒ,
        cₛ = params.cₛ,
    )

    numerical_flux_first_order = model.numerics.flux # should be a function

    dg = DGModel(
        model,
        model.grid.numerical,
        numerical_flux_first_order,
        CentralNumericalFluxSecondOrder(),
        CentralNumericalFluxGradient(),
    )

    return dg
end

coordinates(s::Simulation) = coordinates(simulation.model.grid.numerical)
polynomialorders(s::Simulation) = convention(
    simulation.model.grid.resolution.polynomialorder,
    Val(ndims(simulation.model.grid.domain)),
)

abstract type AbstractTimestepper end

Base.@kwdef struct TimeStepper{S, T} <: AbstractTimestepper
    method::S
    timestep::T
end

"""
calculate_dt(grid, wavespeed = nothing, diffusivity = nothing, viscocity = nothing, cfl = 0.1)
"""
function calculate_dt(
    grid;
    wavespeed = nothing,
    diffusivity = nothing,
    viscocity = nothing,
    cfl = 1.0,
)
    Δx = min_node_distance(grid)
    Δts = []
    if wavespeed != nothing
        push!(Δts, Δx / wavespeed)
    end
    if diffusivity != nothing
        push!(Δts, Δx^2 / diffusivity)
    end
    if viscocity != nothing
        push!(Δts, Δx^2 / viscocity)
    end
    if Δts == []
        @error("Please provide characteristic speed or diffusivities")
        return nothing
    end
    return cfl * minimum(Δts)
end

function calculate_dt(
    grid::DiscretizedDomain;
    wavespeed = nothing,
    diffusivity = nothing,
    viscocity = nothing,
    cfl = 1.0,
)
    return calculate_dt(
        grid.numerical;
        wavespeed = wavespeed,
        diffusivity = diffusivity,
        viscocity = viscocity,
        cfl = cfl,
    )
end
