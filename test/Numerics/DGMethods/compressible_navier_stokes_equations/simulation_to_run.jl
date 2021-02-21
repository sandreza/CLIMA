
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

import ClimateMachine.Mesh.Grids: polynomialorders

using MPI
using LinearAlgebra
using StaticArrays
using Logging, Printf, Dates
# include(pwd() * "/test/Numerics/DGMethods/compressible_navier_stokes_equations/simulation_to_run.jl")
# helper functions
coordinates(s::Simulation) = coordinates(simulation.model.grid.numerical)

polynomialorders(s::Simulation) = convention(simulation.model.grid.resolution.polynomialorder, Val(ndims(simulation.model.grid.domain)))

# initialized on CPU so not problem, but could do kernel launch?
function set_ic!(ϕ, s::Number, x, y, z, p)
    ϕ .= s
    return nothing
end

function set_ic!(ϕ, s::Function, x, y, z, p)
    a_, states, c_ = size(ϕ)
    @inbounds for i in eachindex(x)
        @inbounds for j in 1:states
            ϕʲ = view(ϕ, :, j, :)
            ϕʲ[i] = s(x[i],y[i],z[i],p)[j]
        end
    end
    return nothing
end

function create_callbacks(simulation::Simulation)
    callbacks = simulation.callbacks
    cbvector = [create_callback(simulation, callback) for callback in callbacks]
    return cbvector   
end

function create_callback(simulation::Simulation, callback::Default)
    return nothing
end

config = simulation_to_config(simulation)
dg = config.dg
Q = init_ode_state(dg, Float64(0); init_on_cpu = true)

for s in keys(simulation.initialconditions)
    x, y, z = coordinates(simulation)
    p = simulation.model.parameters
    ic = simulation.initialconditions[s]
    ϕ = getproperty(Q, s)
    set_ic!(ϕ, ic, x, y, z, p)
end

Ns = polynomialorders(simulation)
Nover = (0, 0, 0)
if sum(Nover) > 0
    cutoff = CutoffFilter(dg.grid, Ns .- (Nover .- 1))
    num_state_prognostic = number_states(dg.balance_law, Prognostic())
    Filters.apply!(Q, 1:num_state_prognostic, dg.grid, cutoff)
end

function custom_tendency(tendency, x...; kw...)
    dg(tendency, x...; kw...)
    if sum(Nover) > 0
        cutoff = CutoffFilter(dg.grid, Ns .- (Nover .- 1))
        num_state_prognostic = number_states(dg.balance_law, Prognostic())
        Filters.apply!(tendency, 1:num_state_prognostic, dg.grid, cutoff)
    end
end

Δt = simulation.timestepper.timestep
timestepper = simulation.timestepper.method

odesolver = timestepper(custom_tendency, Q, dt = Δt, t0 = simulation.simulationtime[1])

cbvector = create_callbacks(simulation)

solve!(Q, odesolver; timeend = simulation.simulationtime[2], callbacks = cbvector)

##
BC = (
    ClimateMachine.Ocean.OceanBC(Impenetrable(NoSlip()), Insulating()),
    ClimateMachine.Ocean.OceanBC(
        Impenetrable(KinematicStress(
            (state, aux, t) -> (@SVector [0.01 / state.ρ, -0, -0]),
        )),
        TemperatureFlux((state, aux, t) -> (0.1)),
    ),
)

dg.balance_law.boundary_conditions[2]