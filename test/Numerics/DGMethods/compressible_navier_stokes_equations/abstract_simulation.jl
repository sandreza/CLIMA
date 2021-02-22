import ClimateMachine.Mesh.Grids: polynomialorders

abstract type AbstractSimulation end

Base.@kwdef struct Simulation{𝒜, ℬ, 𝒞, 𝒟, ℰ} <: AbstractSimulation
    model::𝒜
    timestepper::ℬ
    initialconditions::𝒞
    callbacks::𝒟
    simulationtime::ℰ
end

coordinates(s::Simulation) = coordinates(simulation.model.grid.numerical)
polynomialorders(s::Simulation) = convention(simulation.model.grid.resolution.polynomialorder, Val(ndims(simulation.model.grid.domain)))


abstract type AbstractCallback end

struct Default <: AbstractCallback end