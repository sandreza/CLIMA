#!/usr/bin/env julia --project

include("bickley_jet.jl")
ClimateMachine.init()

const FT = Float64

#################
# RUN THE TESTS #
#################

vtkpath =
    abspath(joinpath(ClimateMachine.Settings.output_dir, "vtk_bickley_3D"))

let
    # simulation times
    timeend = FT(200) # s
    dt = FT(0.02) # s
    nout = Int(100)

    # Domain Resolution
    N = 3
    Nˣ = 8
    Nʸ = 8
    Nᶻ = 1

    # Domain size
    Lˣ = 4 * FT(π)  # m
    Lʸ = 4 * FT(π)  # m
    Lᶻ = 1 # m

    params = (; N, Nˣ, Nʸ, Nᶻ, Lˣ, Lʸ, Lᶻ, dt, nout, timeend)

    config = Config(
        "bickley_jet_3D",
        params;
        numerical_flux_first_order = RoeNumericalFlux(),
        Nover = 1,
        periodicity = (true, true, true),
        boundary = ((0, 0), (0, 0), (0, 0)),
        boundary_conditons = (ClimateMachine.Ocean.OceanBC(
            Impenetrable(FreeSlip()),
            Insulating(),
        ),),
    )

    tic = Base.time()

    run_CNSE(config, params; TimeStepper = LSRK54CarpenterKennedy)

    toc = Base.time()
    time = toc - tic
    println(time)
end
