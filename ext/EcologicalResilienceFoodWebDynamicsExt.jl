module EcologicalResilienceFoodWebDynamicsExt

# Datasets of the Domínguez-García et al. (2019) layout from communities of FoodWebDynamics.jl:
# its Jacobian feeds the paper's local metrics, and the simulation experiments run through the
# ResilienceCore community interface that FoodWebDynamics implements.

using EcologicalResilience
using FoodWebDynamics: FoodWebDynamics, Community, biomass, jacobian
using ResilienceCore: PerturbationOptions, global_mortality_tolerance, local_mortality_tolerances, press_response,
                      removal_responses, extinction_tolerance, present_species, state
using LinearAlgebra: Diagonal
using Random: Random, AbstractRNG
import EcologicalResilience: published_pulse_covariance, amplification_metrics, community_record, Jacobian

published_pulse_covariance(c::Community) = Diagonal(c.B[1:FoodWebDynamics.richness(c)] .^ 2)

function amplification_metrics(c::Community; convention::Symbol = :published, dt::Real = 0.1)
    return amplification_metrics(c, jacobian(c.B, c.model), convention, dt)
end
function amplification_metrics(c::Community, J::AbstractMatrix, convention::Symbol, dt::Real)
    convention in (:published, :si) ||
        throw(ArgumentError("convention must be :published or :si, got :$convention"))
    published = convention === :published
    Amax, tmax = max_amplification(J; dt = dt, p = published ? Inf : 2)
    C = published ? published_pulse_covariance(c) : pulse_covariance(biomass(c))
    MAmax, Mtmax = median_max_amplification(J, C; dt = dt)
    return (Amax = Amax, tmax = tmax, MRinf = median_resilience(J, C; dt = dt), MR0 = median_reactivity(J, C),
            MAmax = MAmax, Mtmax = Mtmax)
end

function community_record(rng::AbstractRNG, c::Community; options::PerturbationOptions = PerturbationOptions(),
                          dt::Real = 0.1, convention::Symbol = :published)
    model = c.model
    B = biomass(c)
    N = length(B)
    S = length(state(c))
    J = jacobian(c.B, model)
    amp = amplification_metrics(c, J, convention, dt)
    loc = local_mortality_tolerances(c; options = options)
    rel, l1, _ = press_response(c, fill(1 + options.press, S); options = options)
    pressed(i) = (f = ones(S); f[i] = 1 + options.press; f)
    locpress = [press_response(c, pressed(i); options = options) for i in present_species(c)]
    lrel = first.(locpress)
    ll1 = [x[2] for x in locpress]
    rem = removal_responses(c; options = options)
    surv = [r.survivors / (N - 1) for r in rem]
    rrel = [r.relative for r in rem]
    idx = present_species(c)
    eaten = [any(k -> model.alive[k] && i in model.resources[k], idx) for i in idx]
    p = B ./ sum(B)
    ex = model.exponents
    return (Rinf = resilience(J), R0 = reactivity(J), Is = stochastic_invariability(J), Amax = amp.Amax,
            tmax = amp.tmax, MRinf = amp.MRinf, MR0 = amp.MR0, MAmax = amp.MAmax, Mtmax = amp.Mtmax,
            L1DBRMt_s = press_sensitivity(J), absumS = total_sensitivity(J),
            TM_p = global_mortality_tolerance(c; options = options), minTM_i = minimum(loc),
            av_TM_i = sum(loc) / N, Ext_half = extinction_tolerance(rng, c; options = options), RM_p = rel,
            L1DBRM_p = l1, maxRM_i = lrel[argmax(abs.(lrel))], avRM_i = sum(lrel) / N, avL1DBRM_i = sum(ll1) / N,
            maxSensM_i = maximum(ll1), Robst = sum(surv) / N, DBRobst = sum(rrel) / N,
            L1DBRobst = sum(r.l1 for r in rem) / N, minRobst = minimum(surv), maxDBRost = maximum(rrel),
            minDBRobst = minimum(rrel), maxL1DBRobst = maximum(r.l1 for r in rem), N = N,
            iniN = c.initial_richness, Nb = count(model.basal[idx]), Napex = count(.!eaten),
            Ncan = count(i -> i in model.resources[i], idx), Btot = sum(B), B_IPR = sum(abs2, p),
            a_i = ex.alpha_consumer, a_j = ex.alpha_resource, h_i = ex.eta_consumer, h_j = ex.eta_resource)
end
community_record(c::Community; kwargs...) = community_record(Random.default_rng(), c; kwargs...)

function Jacobian(c::Community)
    idx = present_species(c)
    meta = MatrixMetadata(; taxa = ["sp$i" for i in idx], state_variables = "biomass density",
                          units = "per unit time", derivation = "analytic Jacobian of the bioenergetic model at a polished equilibrium")
    return Jacobian(jacobian(c.B, c.model), biomass(c); meta = meta)
end

end
