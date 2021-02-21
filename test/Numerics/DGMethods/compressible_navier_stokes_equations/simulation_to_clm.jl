# Takes a simulation object and constructs something that works with ClimateMachine
include(pwd() * "/test/Numerics/DGMethods/compressible_navier_stokes_equations/current_template.jl")
include(pwd() * "/test/Numerics/DGMethods/compressible_navier_stokes_equations/three_dimensional/box.jl")

# first construction config

function simulation_to_config(simulation::Simulation; name = "", Nover = 1, mpicomm = MPI.COMM_WORLD, ArrayType = ClimateMachine.array_type())
    numerical_grid = simulation.model.grid.numerical
    dg = simulation_to_model(simulation, simulation.model.balancelaw(), FT = Float64)
    return Config(name, dg, Nover, mpicomm, ArrayType)
end

function simulation_to_model(simulation::Simulation, balancelaw::ThreeDimensionalCompressibleNavierStokesEquations; FT = Float64)
    params = simulation.model.parameters
    if (params.cᶻ == params.cₛ)
        pressure =
            ThreeDimensionalCompressibleNavierStokes.IsotropicPressure{FT}(
                cₛ = params.cₛ,
                ρₒ = params.ρₒ,
            )
    else
        pressure =
            ThreeDimensionalCompressibleNavierStokes.AnisotropicPressure{FT}(
                cₛ = params.cₛ,
                cᶻ = params.cᶻ,
                ρₒ = params.ρₒ,
            )
    end

    Lˣ, Lʸ, Lᶻ = length(simulation.model.grid.domain)

    dissipations = get_dissipation(simulation.model.physics, balancelaw)
    μ = dissipations.μ
    ν = dissipations.ν
    κ = dissipations.κ

    boundary_conditions = get_boundary_conditions(simulation, balancelaw)

    model = ThreeDimensionalCompressibleNavierStokes.CNSE3D(
        (Lˣ, Lʸ, Lᶻ),
        ClimateMachine.Orientations.FlatOrientation(),
        pressure,
        ClimateMachine.Ocean.NonLinearAdvectionTerm(),
        ThreeDimensionalCompressibleNavierStokes.ConstantViscosity{FT}(
            μ = μ,
            ν = ν,
            κ = κ,
        ),
        nothing,
        ThreeDimensionalCompressibleNavierStokes.Buoyancy{FT}(
            α = params.α,
            g = params.g,
        ),
        boundary_conditions,
    )

    numerical_flux_first_order = simulation.model.numerics.flux # should be a function

    dg = DGModel(
        model,
        simulation.model.grid.numerical,
        numerical_flux_first_order,
        CentralNumericalFluxSecondOrder(),
        CentralNumericalFluxGradient(),
    )
    return dg
end


#
function get_dissipation(physics::NamedTuple, ::ThreeDimensionalCompressibleNavierStokesEquations)
    ν = κ = μ = 0
    if haskey(physics, :dissipation)
        if haskey(physics.dissipation, :ρθ)
            κ = physics.dissipation.ρθ.model
        end
        if haskey(physics.dissipation, :ρu)
            ν = physics.dissipation.ρu.model
        end
        if haskey(physics.dissipation, :ρ)
            μ = physics.dissipation.ρ.model
        end
        return (;ν, κ, μ) 
    else
        return (;ν, κ, μ) 
    end
end

get_dissipation(simulation.model.physics, simulation.model.balancelaw())

function get_boundary_conditions(simulation::Simulation, ::ThreeDimensionalCompressibleNavierStokesEquations)
    bcs = simulation.model.boundaryconditions

    westeast   = (check_bc(bcs, :west), check_bc(bcs, :east))
    southnorth = (check_bc(bcs, :south), check_bc(bcs, :north))
    bottomtop  = (check_bc(bcs, :bottom), check_bc(bcs, :top))

    return (westeast, southnorth, bottomtop)
end

function check_bc(bcs, label)
    bctype = ClimateMachine.Ocean.OceanBC
    ρubc = check_bc(bcs, Val(:ρu), label)
    ρθbc = check_bc(bcs, Val(:ρθ), label)
    return bctype(ρubc, ρθbc)
end

function check_bc(bcs, ::Val{:ρθ}, label)
    if haskey(bcs, :ρθ)
        if haskey(bcs[:ρθ], label)
            return bcs[:ρθ][label]
        end
    end
    return Insulating()
end

function check_bc(bcs, ::Val{:ρu}, label)
    if haskey(bcs, :ρu)
        if haskey(bcs[:ρu], label)
            return bcs[:ρu][label]
        end
    end
    return Impenetrable(FreeSlip())
end

import ClimateMachine.Ocean: ocean_init_aux!, ocean_init_state!

function ocean_init_aux!(
    ::ThreeDimensionalCompressibleNavierStokes.CNSE3D,
    aux,
    geom,
)
    @inbounds begin
        aux.x = geom.coord[1]
        aux.y = geom.coord[2]
        aux.z = geom.coord[3]
    end

    return nothing
end

function ocean_init_state!(
    model::ThreeDimensionalCompressibleNavierStokes.CNSE3D,
    state,
    aux,
    localgeo,
    t,
)

    x = aux.x
    y = aux.y
    z = aux.z

    ρ = model.pressure.ρₒ
    state.ρ = ρ
    state.ρu = ρ * @SVector [-0, -0, -0]
    state.ρθ = ρ * 5

    return nothing
end

