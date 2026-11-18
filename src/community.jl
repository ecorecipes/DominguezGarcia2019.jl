# Datasets of the Domínguez-García et al. (2019) layout from communities of FoodWebDynamics.jl:
# the Jacobian feeds ResilienceCore's metrics, and the simulation experiments run through the
# ResilienceCore community interface that FoodWebDynamics implements.

"""
    community_jacobian(c::Community)

The Jacobian of the bioenergetic model at the assembled equilibrium of `c`, as ResilienceCore's
`ContinuousJacobian` with the surviving species' biomasses as its equilibrium (so that the
median metrics' default pulses are proportional to them).
"""
community_jacobian(c::Community) = ContinuousJacobian(jacobian(c), biomass(c))

"""
    published_max_amplification(op::ContinuousJacobian; dt = 0.1)

`(Amax, tmax)` as in Dataset S1: the maximum of the amplification envelope in the ∞-norm
(maximum absolute row sum of `exp(J t)`) on a grid of step `dt`, but `(1.0, 0.0)` whenever the
reactivity in the spectral norm is not positive. In the published data `Amax = 1` exactly when
`R0 ≤ 0`, so the gate uses the 2-norm reactivity `R0` even though the envelope uses the ∞-norm.
ResilienceCore's `max_amplification(op; p = Inf)` gates on the ∞-norm reactivity instead, which
can differ when `R0 ≤ 0` but the ∞-norm logarithmic norm is positive; this helper exists only
to keep the published convention.
"""
function published_max_amplification(op::ContinuousJacobian; dt::Real = 0.1)
    reactivity(op) ≤ 0 && return (1.0, 0.0)
    return max_amplification(op; p = Inf, dt = dt)
end

published_pulse_covariance(c::Community) = Diagonal(c.B[1:FoodWebDynamics.richness(c)] .^ 2)

function amplification_metrics(c::Community; convention::Symbol = :published, dt::Real = 0.1)
    return amplification_metrics(c, community_jacobian(c), convention, dt)
end
function amplification_metrics(c::Community, op::ContinuousJacobian, convention::Symbol, dt::Real)
    convention in (:published, :si) ||
        throw(ArgumentError("convention must be :published or :si, got :$convention"))
    published = convention === :published
    Amax, tmax = published ? published_max_amplification(op; dt = dt) : max_amplification(op; p = 2, dt = dt)
    C = published ? published_pulse_covariance(c) : pulse_covariance(biomass(c))
    MAmax, Mtmax = median_max_amplification(op, C; dt = dt)
    return (Amax = Amax, tmax = tmax, MRinf = median_resilience(op, C; dt = dt), MR0 = median_reactivity(op, C),
            MAmax = MAmax, Mtmax = Mtmax)
end

function community_record(rng::AbstractRNG, c::Community; options::PerturbationOptions = PerturbationOptions(),
                          dt::Real = 0.1, convention::Symbol = :published)
    model = c.model
    B = biomass(c)
    N = length(B)
    S = length(state(c))
    op = community_jacobian(c)
    amp = amplification_metrics(c, op, convention, dt)
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
    return (Rinf = resilience(op), R0 = reactivity(op), Is = stochastic_invariability(op), Amax = amp.Amax,
            tmax = amp.tmax, MRinf = amp.MRinf, MR0 = amp.MR0, MAmax = amp.MAmax, Mtmax = amp.Mtmax,
            L1DBRMt_s = press_sensitivity(op), absumS = total_sensitivity(op),
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
