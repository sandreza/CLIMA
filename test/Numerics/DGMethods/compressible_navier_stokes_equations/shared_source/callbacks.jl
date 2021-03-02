abstract type AbstractCallback end

struct Default <: AbstractCallback end
struct Info <: AbstractCallback end
struct StateCheck{T} <: AbstractCallback
    number_of_checks::T
end

@Base.kwdef struct JLD2State{T, V, B} <: AbstractCallback
    iteration::T
    filepath::V
    overwrite::B = true
end

function create_callback(output::JLD2State, simulation::Simulation, odesolver)
    # Initialize output
    output.overwrite && isfile(output.filepath) && rm(output.filepath; force = output.overwrite)

    Q = simulation.state
    mpicomm = MPI.COMM_WORLD
    iteration = output.iteration

    steps = ClimateMachine.ODESolvers.getsteps(odesolver)
    time = ClimateMachine.ODESolvers.gettime(odesolver)

    file = jldopen(output.filepath, "a+")
    JLD2.Group(file, "state")
    JLD2.Group(file, "time")
    file["state"][string(steps)] = Array(Q)
    file["time"][string(steps)] = time
    close(file)


    jldcallback = ClimateMachine.GenericCallbacks.EveryXSimulationSteps(
        iteration
    ) do (s = false)
        steps = ClimateMachine.ODESolvers.getsteps(odesolver)
        time = ClimateMachine.ODESolvers.gettime(odesolver)
        @info steps, time
        file = jldopen(output.filepath, "a+")
        file["state"][string(steps)] = Array(Q)
        file["time"][string(steps)] = time
        close(file)
        return nothing
    end
end

function create_callbacks(simulation::Simulation, odesolver)
    callbacks = simulation.callbacks

    if isempty(callbacks)
        return ()
    else
        cbvector = [
            create_callback(callback, simulation, odesolver)
            for callback in callbacks
        ]
        return tuple(cbvector...)
    end
end

function create_callback(::Default, simulation::Simulation, odesolver)
    cb_info = create_callback(Info(), simulation, odesolver)
    cb_state_check = create_callback(StateCheck(10), simulation, odesolver)

    return (cb_info, cb_state_check)
end

function create_callback(::Info, simulation::Simulation, odesolver)
    Q = simulation.state
    timeend = simulation.time.finish
    mpicomm = MPI.COMM_WORLD

    starttime = Ref(now())
    cbinfo = ClimateMachine.GenericCallbacks.EveryXWallTimeSeconds(
        60,
        mpicomm,
    ) do (s = false)
        if s
            starttime[] = now()
        else
            energy = norm(Q)
            @info @sprintf(
                """Update
                simtime = %8.2f / %8.2f
                runtime = %s
                norm(Q) = %.16e""",
                ClimateMachine.ODESolvers.gettime(odesolver),
                timeend,
                Dates.format(
                    convert(Dates.DateTime, Dates.now() - starttime[]),
                    Dates.dateformat"HH:MM:SS",
                ),
                energy
            )

            if isnan(energy)
                error("NaNs")
            end
        end
    end

    return cbinfo
end

function create_callback(callback::StateCheck, simulation::Simulation, _...)
    sim_length = simulation.time.finish - simulation.time.start
    timestep = simulation.timestepper.timestep
    nChecks = callback.number_of_checks


    nt_freq = floor(Int, sim_length / timestep / nChecks)

    cbcs_dg = ClimateMachine.StateCheck.sccreate(
        [(simulation.state, "state")],
        nt_freq,
    )

    return cbcs_dg
end
