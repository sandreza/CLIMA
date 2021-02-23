abstract type AbstractCallback end

struct Default <: AbstractCallback end
struct Info <: AbstractCallback end

function create_callbacks(simulation::Simulation)
    callbacks = simulation.callbacks
    cbvector = [create_callback(simulation, callback) for callback in callbacks]
    return cbvector
end

function create_callback(simulation::Simulation, callback::Default)
    return nothing
end

function create_callback(simulation::Simulation, callback::Info)
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
                ODESolvers.gettime(odesolver),
                timespan.timeend,
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
