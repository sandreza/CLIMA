include("../CNSE.jl")
include("ThreeDimensionalCompressibleNavierStokesEquations.jl")
# Config("", (Nˣ = 1, Nʸ=1, Nᶻ = 1, N = 1), (Lˣ = 1, Lʸ = 1, Lᶻ = 1), (cₛ=1, cᶻ=1, ρₒ=1, μ=0, ν=0, κ=0, α=0, g=0))
function Config(
    name,
    resolution,
    domain,
    params;
    numerical_flux_first_order = RusanovNumericalFlux(),
    Nover = 0,
    periodicity = (true, true, true),
    boundary = ((0, 0), (0, 0), (0, 0)),
    boundary_conditions = (),
)
    mpicomm = MPI.COMM_WORLD
    ArrayType = ClimateMachine.array_type()

    xrange =
        range(-domain.Lˣ / 2; length = resolution.Nˣ + 1, stop = domain.Lˣ / 2)
    yrange =
        range(-domain.Lʸ / 2; length = resolution.Nʸ + 1, stop = domain.Lʸ / 2)
    zrange =
        range(-domain.Lᶻ / 2; length = resolution.Nᶻ + 1, stop = domain.Lᶻ / 2)

    brickrange = (xrange, yrange, zrange)

    topl = BrickTopology(
        mpicomm,
        brickrange,
        periodicity = periodicity,
        boundary = boundary,
    )

    grid = DiscontinuousSpectralElementGrid(
        topl,
        FloatType = FT,
        DeviceArray = ArrayType,
        polynomialorder = resolution.N + Nover,
    )

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

    model = ThreeDimensionalCompressibleNavierStokes.CNSE3D(
        (domain.Lˣ, domain.Lʸ, domain.Lᶻ),
        ClimateMachine.Orientations.FlatOrientation(),
        pressure,
        ClimateMachine.Ocean.NonLinearAdvectionTerm(),
        ThreeDimensionalCompressibleNavierStokes.ConstantViscosity{FT}(
            μ = params.μ,
            ν = params.ν,
            κ = params.κ,
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

    return Config(name, dg, Nover, mpicomm, ArrayType)
end

import ClimateMachine.Ocean: ocean_init_aux!

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
