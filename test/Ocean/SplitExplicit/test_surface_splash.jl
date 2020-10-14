#!/usr/bin/env julia --project

include("hydrostatic_spindown.jl")
ClimateMachine.init()

const FT = Float64

#################
# RUN THE TESTS #
#################

import ClimateMachine.Ocean: ocean_init_state!

function ocean_init_state!(
    m::Union{HydrostaticBoussinesqModel, OceanModel},
    p::SimpleBox,
    Q,
    A,
    coords,
    t,
)
    @inbounds x = coords[1]
    @inbounds y = coords[2]

    u = -0
    v = -0
    η = exp(-(x^2 + y^2) / (p.Lˣ^2 + p.Lʸ^2))

    Q.u = @SVector[u, v]
    Q.η = η
    Q.θ = -0

    return nothing
end

function ocean_init_state!(
    m::Union{ShallowWaterModel, BarotropicModel},
    p::SimpleBox,
    Q,
    A,
    coords,
    t,
)
    @inbounds x = coords[1]
    @inbounds y = coords[2]

    U = -0
    V = -0
    η = exp(-(x^2 + y^2) / (p.Lˣ^2 + p.Lʸ^2))

    Q.U = @SVector[U, V]
    Q.η = η

    return nothing
end

@testset "$(@__FILE__)" begin

    # simulation time
    timeend = FT(24 * 3600) # s
    tout = FT(1.5 * 3600) # s
    timespan = (tout, timeend)

    # DG polynomial order
    N = Int(4)

    # Domain resolution
    Nˣ = Int(5)
    Nʸ = Int(5)
    Nᶻ = Int(8)
    resolution = (N, Nˣ, Nʸ, Nᶻ)

    # Domain size
    Lˣ = 1e6  # m
    Lʸ = 1e6  # m
    H = 400  # m
    dimensions = (Lˣ, Lʸ, H)

    problem = SimpleBox{FT}(
        dimensions...;
        BC = (
            OceanBC(Impenetrable(FreeSlip()), Insulating()),
            OceanBC(Penetrable(FreeSlip()), Insulating()),
        ),
        rotation = Fixed(),
    )

    config = SplitConfig(
        "surface_splash",
        resolution,
        dimensions,
        Coupled(),
        problem;
        solver = SplitExplicitSolver,
    )

    run_split_explicit(
        config,
        timespan;
        dt_fast = 300, # seconds
        dt_slow = 300, # 90 * 60, # seconds
        analytic_solution = true,
    )
end