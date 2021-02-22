#######
# useful concepts for dispatch
#######

"""
Advection terms

right now really only non-linear or ::Nothing
"""
abstract type AdvectionTerm end
struct NonLinearAdvectionTerm <: AdvectionTerm end

"""
Turbulence Closures

ways to handle drag and diffusion and such
"""
abstract type TurbulenceClosure end

struct LinearDrag{T} <: TurbulenceClosure
    λ::T
end

struct ConstantViscosity{T} <: TurbulenceClosure
    μ::T
    ν::T
    κ::T
    function ConstantViscosity{T}(;
        μ = T(1e-6),   # m²/s
        ν = T(1e-6),   # m²/s
        κ = T(1e-6),   # m²/s
    ) where {T <: AbstractFloat}
        return new{T}(μ, ν, κ)
    end
end

"""
Forcings

ways to add body terms and sources
"""
abstract type Forcing end
abstract type CoriolisForce <: Forcing end

struct fPlaneCoriolis{T} <: CoriolisForce
    fₒ::T
    β::T
    function fPlaneCoriolis{T}(;
        fₒ = T(1e-4), # Hz
        β = T(1e-11), # Hz/m
    ) where {T <: AbstractFloat}
        return new{T}(fₒ, β)
    end
end

struct WindStress{T} <: Forcing
    τₒ::T
    function WindStress{T}(; τₒ = T(1e-4)) where {T <: AbstractFloat}
        return new{T}(τₒ)
    end
end

struct Buoyancy{T} <: Forcing
    α::T # 1/K
    g::T # m/s²
    function Buoyancy{T}(; α = T(2e-4), g = T(10)) where {T <: AbstractFloat}
        return new{T}(α, g)
    end
end


