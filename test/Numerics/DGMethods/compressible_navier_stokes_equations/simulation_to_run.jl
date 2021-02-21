
using Test
using JLD2
using ClimateMachine
ClimateMachine.init()
using ClimateMachine.ODESolvers
using ClimateMachine.Mesh.Filters
using ClimateMachine.VariableTemplates
using ClimateMachine.Mesh.Grids: polynomialorders
using ClimateMachine.Ocean

using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.BalanceLaws:
    vars_state, Prognostic, Auxiliary, number_states
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.MPIStateArrays
using ClimateMachine.VTK

using MPI
using LinearAlgebra
using StaticArrays
using Logging, Printf, Dates

function run_CNSE(
    config,
    resolution,
    timespan;
    TimeStepper = LSRK54CarpenterKennedy,
    refDat = (),
)
    dg = config.dg
    Q = init_ode_state(dg, FT(0); init_on_cpu = true)

    if config.Nover > 0
        cutoff = CutoffFilter(dg.grid, resolution.N + config.Nover)
        num_state_prognostic = number_states(dg.balance_law, Prognostic())
        Filters.apply!(Q, 1:num_state_prognostic, dg.grid, cutoff)
    end

    function custom_tendency(tendency, x...; kw...)
        dg(tendency, x...; kw...)
        if config.Nover > 0
            cutoff = CutoffFilter(dg.grid, resolution.N + config.Nover)
            num_state_prognostic = number_states(dg.balance_law, Prognostic())
            Filters.apply!(tendency, 1:num_state_prognostic, dg.grid, cutoff)
        end
    end

    println("time step is " * string(timespan.dt))
    println("time end is " * string(timespan.timeend))
    odesolver = TimeStepper(custom_tendency, Q, dt = timespan.dt, t0 = 0)

    vtkstep = 0
    cbvector = make_callbacks(
        vtkpath,
        vtkstep,
        timespan,
        config.mpicomm,
        odesolver,
        dg,
        dg.balance_law,
        Q,
        filename = config.name,
    )

    eng0 = norm(Q)
    @info @sprintf """Starting
    norm(Q₀) = %.16e
    ArrayType = %s""" eng0 config.ArrayType

    solve!(Q, odesolver; timeend = timespan.timeend, callbacks = cbvector)

    ## Check results against reference
    ClimateMachine.StateCheck.scprintref(cbvector[end])
    if length(refDat) > 0
        @test ClimateMachine.StateCheck.scdocheck(cbvector[end], refDat)
    end

    return nothing
end