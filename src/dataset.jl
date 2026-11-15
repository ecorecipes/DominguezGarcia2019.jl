# One row per community in the layout of the authors' raw Dataset S1/S2 (GitHub deposit), and
# the transformation to the paper's 27 "larger is more stable" metrics performed by
# `prepare_dataframe` in the authors' Functions.py.

"""
    RAW_COLUMNS

The 39 data columns of the raw Dataset S1/S2 (after the `network` identifier), in file order:
28 metric columns, then `N`, `iniN`, `Nb`, `Napex`, `Ncan`, `Btot`, `B_IPR` and the
community's exponents `a_i`, `a_j`, `h_i`, `h_j`.
"""
const RAW_COLUMNS = [:Rinf, :R0, :Is, :Amax, :tmax, :MRinf, :MR0, :MAmax, :Mtmax, :L1DBRMt_s, :absumS, :TM_p,
                     :minTM_i, :av_TM_i, :Ext_half, :RM_p, :L1DBRM_p, :maxRM_i, :avRM_i, :avL1DBRM_i,
                     :maxSensM_i, :Robst, :DBRobst, :L1DBRobst, :minRobst, :maxDBRost, :minDBRobst,
                     :maxL1DBRobst, :N, :iniN, :Nb, :Napex, :Ncan, :Btot, :B_IPR, :a_i, :a_j, :h_i, :h_j]

"""
    METRICS

The labels of the 27 stability metrics, in the order of `measures` in the authors' code.
"""
const METRICS = ["R0", "Is", "Amax", "tmax", "Rinf", "MR0", "MAmax", "Mtmax", "MRinf", "S", "<s_ij>", "TM",
                 "<TM>", "minTM", "<CE>", "TE", "maxCE", "RM", "<RM>", "maxRM", "SM", "<SM>", "maxSM", "<RE>",
                 "maxRE", "<SE>", "maxSE"]

# Publisher (PNAS) column labels of the same 27 metrics, in the same order.
const PUBLISHER_LABELS = ["\$R_{0}\$", "\$Is\$", "\$A_{max}\$", "\$t_{max}\$", "\$R_{inf}\$", "\$MR_{0}\$",
                          "\$MA_{max}\$", "\$Mt_{max}\$", "\$MR_{inf}\$", "\$S\$", "\$<s_{ij}>\$", "\$TM^G\$",
                          "\$<TM^L>\$", "\$TM^L_{min}\$", "\$<CE>\$", "\$<TE>\$", "\$CE_{max}\$", "\$RM^G\$",
                          "\$<RM^L>\$", "\$RM^L_{max}\$", "\$SM^G\$", "\$<SM^L>\$", "\$SM^L_{max}\$",
                          "\$<RE>\$", "\$RE_{max}\$", "\$<SE>\$", "\$SE_{max}\$"]

"""
    published_pulse_covariance(c::Community)

The pulse covariance behind the published median metrics (`MR0`, `MRinf`, `MAmax`, `Mtmax`) of
Dataset S1: `Diagonal(w.^2)` with `w = c.B[1:N]`. That is the first `N` entries of the biomass
vector over all initial species (extinct ones at zero, in niche order; see
[`niche_model`](@ref)), set against the `N` surviving species. The weights agree with
[`pulse_covariance`](@ref)`(biomass(c))` only when no species went extinct during assembly.
The misalignment reproduces the published data: in communities that kept every species the
published median metrics match the aligned weights, with `MR0 ≈ −0.76` and no amplification.
After extinctions, `MR0` drifts towards 0 and about a fifth of communities amplify, as in the
published data; aligned weights give neither.
"""
published_pulse_covariance(c::Community) = Diagonal(c.B[1:richness(c)] .^ 2)

"""
    amplification_metrics(c::Community; convention = :published, dt = 0.1)

The six raw columns of [`community_record`](@ref) that depend on the convention, as a
`NamedTuple` `(Amax, tmax, MRinf, MR0, MAmax, Mtmax)`.

- `:published` reproduces Dataset S1. `Amax` and `tmax` come from the envelope in the ∞-norm
  ([`max_amplification`](@ref) with `p = Inf`), and the median metrics use
  [`published_pulse_covariance`](@ref).
- `:si` follows the definitions of the SI: the spectral norm, and pulses proportional to the
  surviving species' own biomasses ([`pulse_covariance`](@ref)).
"""
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

"""
    community_record([rng], c::Community; options = PerturbationOptions(), dt = 0.1,
                     convention = :published)

All raw columns of Dataset S1 ([`RAW_COLUMNS`](@ref)) for one assembled community, as a
`NamedTuple`. This runs every perturbation experiment, so it is the expensive step.
`convention` selects how `Amax`, `tmax` and the four median metrics are computed (see
[`amplification_metrics`](@ref)). `:published` reproduces Dataset S1, and `:si` follows the
SI's definitions.

The global press (`RM_p`, `L1DBRM_p`) raises every species' mortality by `options.press`, as
the local presses do. Ranked within each richness, the published `RM_p` agrees with the sum of
the local responses (Spearman ρ ≈ 0.92), as this definition does, but its values are about 30
times larger. No stronger press tried (mortality ×2–×11, plant carrying capacity, all loss
terms) reproduces both the size and the ranking. The scale does not affect the paper's
rank-based analyses.
"""
function community_record(rng::AbstractRNG, c::Community; options::PerturbationOptions = PerturbationOptions(),
                          dt::Real = 0.1, convention::Symbol = :published)
    model = c.model
    B = biomass(c)
    N = length(B)
    J = jacobian(c.B, model)
    amp = amplification_metrics(c, J, convention, dt)
    loc = local_mortality_tolerances(c; options = options)
    rel, l1, _ = press_response(c, (m, _) -> set_mortality!(m, 1 + options.press); options = options)
    locpress = [press_response(c, (m, _) -> (m.mortality[i] = 1 + options.press); options = options)
                for i in findall(model.alive)]
    lrel = first.(locpress)
    ll1 = [x[2] for x in locpress]
    rem = removal_responses(c; options = options)
    surv = [r.survivors / (N - 1) for r in rem]
    rrel = [r.relative for r in rem]
    idx = findall(model.alive)
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

"""
    prepare_metrics(raw::DataFrame)

The 27 stability metrics ([`METRICS`](@ref)) from raw Dataset-S1-style columns, transformed as
the authors' `prepare_dataframe` does so that larger values mean more stability: sign flips,
absolute values, normalisation of `absumS` by N² and of `Ext_half` by N, and the secondary
extinction counts `N(1 − Robst)`. Returns a `DataFrame` with `network` (if present), the 27
metrics and `N`. On the authors' data this equals the publisher's Dataset S1/S2 column by column.
"""
function prepare_metrics(raw::DataFrame)
    N = Float64.(raw.N)
    out = DataFrame()
    hasproperty(raw, :network) && (out.network = raw.network)
    cols = (R0 = -raw.R0, Is = raw.Is, Amax = -raw.Amax, tmax = -raw.tmax, Rinf = raw.Rinf, MR0 = -raw.MR0,
            MAmax = -raw.MAmax, Mtmax = -raw.Mtmax, MRinf = raw.MRinf, S = -raw.L1DBRMt_s,
            s_ij = -raw.absumS ./ N .^ 2, TM = raw.TM_p, TM_mean = raw.av_TM_i, minTM = raw.minTM_i,
            CE = -N .* (1 .- raw.Robst), TE = raw.Ext_half ./ N, maxCE = -N .* (1 .- raw.minRobst),
            RM = -abs.(raw.RM_p), RM_mean = -abs.(raw.avRM_i), maxRM = -abs.(raw.maxRM_i), SM = -raw.L1DBRM_p,
            SM_mean = -raw.avL1DBRM_i, maxSM = -raw.maxSensM_i, RE = -abs.(raw.DBRobst),
            maxRE = -max.(abs.(raw.maxDBRost), abs.(raw.minDBRobst)), SE = -raw.L1DBRobst,
            maxSE = -raw.maxL1DBRobst)
    for (label, v) in zip(METRICS, values(cols))
        out[!, label] = Float64.(v)
    end
    out.N = raw.N
    return out
end

"""
    read_metrics(path)

Read a dataset of stability metrics and return it in the form of [`prepare_metrics`](@ref):
either a raw Dataset S1/S2 from the authors' GitHub deposit or one written by
[`community_record`](@ref) (both have the column `L1DBRMt_s`), or the publisher's Dataset S1/S2
(LaTeX column labels, already transformed).
"""
function read_metrics(path::AbstractString)
    df = CSV.read(path, DataFrame)
    hasproperty(df, :L1DBRMt_s) && return prepare_metrics(df)
    all(l -> l in names(df), PUBLISHER_LABELS) ||
        throw(ArgumentError("$path is neither a raw nor a publisher-format metrics dataset"))
    out = DataFrame(network = df.network)
    for (label, pub) in zip(METRICS, PUBLISHER_LABELS)
        out[!, label] = Float64.(df[!, pub])
    end
    out.N = df.N
    return out
end
