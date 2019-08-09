#### Moisture component in atmosphere model
abstract type MoistureModel end

vars_state(::MoistureModel, T) = Tuple{}
vars_gradient(::MoistureModel, T) = Tuple{}
vars_diffusive(::MoistureModel, T) = Tuple{}
vars_aux(::MoistureModel, T) = Tuple{}

function update_aux!(::MoistureModel, state::Vars, diffusive::Vars, aux::Vars, t::Real)
end
function diffusive!(::MoistureModel, diffusive, ∇transform, state, aux, t, ν)
end
function flux_diffusive!(::MoistureModel, flux::Grad, state::Vars, diffusive::Vars, aux::Vars, t::Real)
end
function gradvariables!(::MoistureModel, transform::Vars, state::Vars, aux::Vars, t::Real)
end

function internal_energy(m::MoistureModel, state::Vars, aux::Vars)
  T = eltype(state)
  ρinv = 1 / state.ρ
  ρe_kin = ρinv*sum(abs2, state.ρu)/2
  ρe_pot = state.ρ * grav * aux.coord.z
  ρe_int = state.ρe - ρe_kin - ρe_pot
  e_int = ρinv*ρe_int
  return e_int
end

temperature(m::MoistureModel, state::Vars, aux::Vars) = air_temperature(thermo_state(m, state, aux))
pressure(m::MoistureModel, state::Vars, aux::Vars) = air_pressure(thermo_state(m, state, aux))
soundspeed(m::MoistureModel, state::Vars, aux::Vars) = soundspeed_air(thermo_state(m, state, aux))

"""
    DryModel

Assumes the moisture components is in the dry limit.
"""
struct DryModel <: MoistureModel
end

vars_aux(::DryModel,T) = NamedTuple{(:e_int, :temperature), Tuple{T,T}}
function update_aux!(m::DryModel, state::Vars, diffusive::Vars, aux::Vars, t::Real)
  aux.moisture.e_int = internal_energy(m, state, aux)
  TS = PhaseDry(aux.moisture.e_int, state.ρ)
  aux.moisture.temperature = air_temperature(TS)
  nothing
end

thermo_state(::DryModel, state::Vars, aux::Vars) = PhaseDry(aux.moisture.e_int, state.ρ)

"""
    EquilMoist

Assumes the moisture components are computed via thermodynamic equilibrium.
"""
struct EquilMoist <: MoistureModel
end
vars_state(::EquilMoist,T) = NamedTuple{(:ρq_tot,),Tuple{T}}
vars_gradient(::EquilMoist,T) = NamedTuple{(:q_tot, :total_enthalpy),Tuple{T,T,T,T}}
vars_diffusive(::EquilMoist,T) = NamedTuple{(:ρd_q_tot, :ρ_SGS_enthalpyflux), Tuple{SVector{3,T},SVector{3,T}}}
vars_aux(::EquilMoist,T) = NamedTuple{(:e_int, :temperature), Tuple{T,T}}

function update_aux!(m::EquilMoist, state::Vars, diffusive::Vars, aux::Vars, t::Real)
  aux.moisture.e_int = internal_energy(m, state, aux)
  TS = PhaseEquil(aux.moisture.e_int, get_phase_partition(m, state).tot, state.ρ)
  aux.moisture.temperature = air_temperature(TS)
  nothing
end


get_phase_partition(::EquilMoist, state::Vars) = PhasePartition(state.ρq_tot/state.ρ)
thermo_state(::EquilMoist, state::Vars, aux::Vars) = PhaseEquil(aux.moisture.e_int, state.ρq_tot/state.ρ, state.ρ, aux.moisture.temperature)

function gradvariables!(m::EquilMoist, transform::Vars, state::Vars, aux::Vars, t::Real)
  invρ = state.ρ
  transform.moisture.q_tot = state.ρq_tot * invρ

  phase = thermo_state(m, state, aux)
  q = PhasePartition(phase)
  R_m = gas_constant_air(q)
  T = aux.moisture.temperature
  e_tot = state.ρe * invρ
  transform.moisture.total_enthalpy = e_tot + R_m*T
end


function diffusive!(m::EquilMoist, diffusive::Vars, ∇transform::Grad, state::Vars, aux::Vars, t::Real, ν::Union{Real,AbstractMatrix})
  # turbulent Prandtl number
  diag_ν = ν isa Real ? ν : diag(ν) # either a scalar or vector
  D_T = diag_ν / Prandtl_t

  # diffusive flux of q_tot
  diffusive.moisture.ρd_q_tot = state.ρ * (-D_T) .* ∇transform.moisture.q_tot

  # diffusive flux of total energy
  diffusive.moisture.ρ_SGS_enthalpyflux = state.ρ * (-D_T) .* ∇transform.transform.moisture.total_enthalpy
end

function flux_diffusive!(m::EquilMoist, flux::Grad, state::Vars, diffusive::Vars, aux::Vars, t::Real)
  flux.ρ += diffusive.moisture.ρd_q_tot
  flux.ρu += diffusive.moisture.ρd_q_tot .* u'
  flux.ρe += diffusive.moisture.ρJ_ρD

  flux.moisture.ρq_tot = diffusive.moisture.ρd_q_tot
end
