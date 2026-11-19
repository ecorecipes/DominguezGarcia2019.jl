"""
    DominguezGarcia2019

A reproduction of Domínguez-García, Dakos & Kéfi (2019), "Unveiling dimensions of stability in
complex ecological networks", PNAS 116:25714 (doi:10.1073/pnas.1904470116).

This is a case study: it holds only what is specific to the paper.

- [`community_record`](@ref) computes the 39 raw columns of the paper's Dataset S1 for an
  assembled food web, with the conventions the published data actually used
  ([`amplification_metrics`](@ref), [`published_pulse_covariance`](@ref),
  [`published_max_amplification`](@ref)).
- [`prepare_metrics`](@ref) and [`read_metrics`](@ref) turn raw columns into the paper's 27
  "larger is more stable" metrics.
- The analysis functions ([`size_correlations`](@ref), [`modular_groups`](@ref),
  [`explained_variance`](@ref), …) port the authors' Python code.

The metrics themselves come from ResilienceCore.jl, on a `ContinuousJacobian` (see
[`community_jacobian`](@ref)), and are re-exported. The food webs, the bioenergetic model and
community assembly come from FoodWebDynamics.jl. Interaction matrices, structural stability
and sign patterns are in CommunityMatrices.jl.
"""
module DominguezGarcia2019

using CSV: CSV
using Clustering: Hclust, hclust
using DataFrames: DataFrame, groupby, rename!
using LinearAlgebra: Diagonal, Symmetric, diag, diagind, eigvals, norm, tr
using Statistics: cor, cov, mean, std
using StatsBase: StatsBase
using Random: Random, AbstractRNG
using FoodWebDynamics: FoodWebDynamics, Community, biomass, jacobian
using ResilienceCore: ResilienceCore, ContinuousJacobian, resilience, reactivity, amplification_envelope,
                      max_amplification, pulse_covariance, median_reactivity, median_displacement,
                      median_max_amplification, median_resilience, stochastic_invariability, sensitivity_matrix,
                      press_sensitivity, total_sensitivity, PerturbationOptions, global_mortality_tolerance,
                      local_mortality_tolerances, press_response, removal_responses, extinction_tolerance,
                      present_species, state, global_press_response,
                      local_press_responses

# metrics and experiments of ResilienceCore, re-exported for the vignettes and scripts
export ContinuousJacobian, resilience, reactivity, amplification_envelope, max_amplification, pulse_covariance,
       median_reactivity, median_displacement, median_max_amplification, median_resilience,
       stochastic_invariability, sensitivity_matrix, press_sensitivity, total_sensitivity
export PerturbationOptions, global_mortality_tolerance, local_mortality_tolerances, press_response,
       removal_responses, extinction_tolerance
# the paper's datasets and conventions
export RAW_COLUMNS, METRICS, community_jacobian, published_max_amplification, published_pulse_covariance,
       amplification_metrics, community_record, prepare_metrics, read_metrics
# the paper's analysis
export SizeCorrelations, size_correlations, class_average, SIZE_CLASSES, metric_network, modularity,
       modular_groups, cophenetic, metric_dendrogram, explained_variance, ellipsoid_volume, correlation_types

include("dataset.jl")
include("community.jl")
include("analysis.jl")

end
