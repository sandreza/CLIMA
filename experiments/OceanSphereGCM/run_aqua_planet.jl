#!/usr/bin/env julia --project

include("aqua_planet.jl")
ClimateMachine.init()

# Float type
const FT = Float64

let
    # simulation time
    timestart = FT(0)      # s
    timestep = FT(1)     # s
    timeend = FT(1) # s
    timespan = (timestart, timeend)

    # DG polynomial order
    N = Int(4)

    # Domain resolution
    Nʰ = Int(12)
    Nᶻ = Int(6)
    resolution = (N, Nʰ, Nᶻ)

    # Domain size
    domain_height = FT(1000)   # m

    BC = (
        OceanBC(Impenetrable(NoSlip()), Insulating()),
        OceanBC(Penetrable(FreeSlip()), Insulating()),
    )

    config = config_aqua_planet(
        "aqua_planet",
        resolution,
        domain_height,
        SimpleSphere,
        BC,
    )

    run_aqua_planet(config, timespan, timestep)
end
