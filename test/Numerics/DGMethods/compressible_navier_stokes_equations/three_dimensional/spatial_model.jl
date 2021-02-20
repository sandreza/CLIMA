abstract type AbstractModel end

struct SpatialModel{𝒜,ℬ,𝒞,𝒟,ℰ,ℱ,𝒢} <: AbstractModel
    balancelaw::𝒜
    domain::ℬ
    boundaryconditions::𝒞
    fluxes::𝒟
    physics::ℰ
    dissipation::ℱ
    parameters::𝒢
end