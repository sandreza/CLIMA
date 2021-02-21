# Takes a simulation object and constructs something that works with ClimateMachine
include(pwd() * "/test/Numerics/DGMethods/compressible_navier_stokes_equations/current_template.jl")

# this does not seem like the right idea
ThreeDimensionalCompressibleNavierStokesEquations() = ThreeDimensionalCompressibleNavierStokesEquations(1,1,1,1,1,1,1,1,1)
# first construction config

function simulation_to_config(simulation::Simulation; name = "", Nover = 1, mpicomm = MPI.COMM_WORLD, ArrayType = ClimateMachine.array_type())
    println("hello")
    numerical_grid = simulation.model.grid.numerical
    dg = simulation_to_model(simulation, simulation.model.balancelaw(), FT = eltype(grid))
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

    model = ThreeDimensionalCompressibleNavierStokesEquations(
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

    dg = DGModel(
        model,
        grid,
        numerical_flux_first_order,
        CentralNumericalFluxSecondOrder(),
        CentralNumericalFluxGradient(),
    )
    return dg

end


##
function get_dissipation(physics::NamedTuple, ::ThreeDimensionalCompressibleNavierStokesEquations)
    ν = κ = μ = 0
    if :dissipation in keys(simulation.model.physics)
        if :ρθ in keys(simulation.model.physics.dissipation)
            κ = simulation.model.physics.dissipation.ρθ.model
        end
        if :ρu in keys(simulation.model.physics.dissipation)
            ν = simulation.model.physics.dissipation.ρu.model
        end
        if :ρ in keys(simulation.model.physics.dissipation)
            μ = simulation.model.physics.dissipation.ρ.model
        end
        return (;ν, κ, μ) 
    else
        return (;ν, κ, μ) 
    end
end

get_dissipation(simulation.model.physics, simulation.model.balancelaw())

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

