abstract type AbstractModel end

Base.@kwdef struct SpatialModel{𝒜,ℬ,𝒞,𝒟,ℰ,ℱ} <: AbstractModel
    balancelaw::𝒜
    physics::ℬ
    numerics::𝒞
    grid::𝒟
    boundaryconditions::ℰ
    parameters::ℱ
end