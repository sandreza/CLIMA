abstract type AbstractSimulation end

Base.@kwdef struct Simulation{𝒜, ℬ, 𝒞, 𝒟, ℰ} <: AbstractSimulation
    model::𝒜
    timestepper::ℬ
    initialconditions::𝒞
    callbacks::𝒟
    simulationtime::ℰ
end


abstract type AbstractCallback end

struct Default <: AbstractCallback end