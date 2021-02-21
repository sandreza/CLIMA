abstract type AbstractTimestepper end

Base.@kwdef struct TimeStepper{S, T} <: AbstractTimestepper
    method::S
    timestep::T
end

"""
calculate_dt(grid, wavespeed = nothing, diffusivity = nothing, viscocity = nothing, cfl = 0.1)
"""
function calculate_dt(grid; wavespeed = nothing, diffusivity = nothing, viscocity = nothing, cfl = 1.0)
    Δx = min_node_distance(grid)
    Δts = []
    if wavespeed != nothing
        push!(Δts, Δx/wavespeed)
    end
    if diffusivity != nothing
        push!(Δts, Δx^2/diffusivity)
    end
    if viscocity != nothing
        push!(Δts, Δx^2/viscocity)
    end
    if Δts == []
        @error("Please provide characteristic speed or diffusivities")
        return nothing
    end
    return cfl * minimum(Δts)
end

calculate_dt(grid::DiscretizedDomain; wavespeed = nothing, diffusivity = nothing, viscocity = nothing, cfl = 1.0) = calculate_dt(grid.numerical; wavespeed = wavespeed, diffusivity = diffusivity, viscocity = viscocity, cfl = cfl)
