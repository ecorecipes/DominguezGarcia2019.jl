"""
    EcologicalResilience

Resilience and stability metrics for ecological models and data, starting with an end-to-end
reproduction of Domínguez-García, Dakos & Kéfi (2019), "Unveiling dimensions of stability in
complex ecological networks", PNAS 116:25714 (doi:10.1073/pnas.1904470116): niche-model food
webs, an allometric bioenergetic model, 27 stability metrics, and the correlation, modularity,
clustering and variance analyses of the paper. Plans and verification notes: the OKF wiki in
`../../wiki`.
"""
module EcologicalResilience

using CSV: CSV
using Clustering: Hclust, hclust
using DataFrames: DataFrame, groupby, rename!
using LinearAlgebra: I, Diagonal, Symmetric, diag, diagind, dot, eigmax, eigvals, issuccess, kron, lu, mul!, norm,
                     opnorm, svdvals, sylvester, tr
using Statistics: cor, cov, mean, std
using StatsBase: StatsBase
using OrdinaryDiffEqTsit5: Tsit5
using Random: Random, AbstractRNG
using SciMLBase: SciMLBase, ODEFunction, ODEProblem, DiscreteCallback
using ResilienceCore: ResilienceCore, PerturbationOptions, global_mortality_tolerance, local_mortality_tolerances,
                      press_response, removal_responses, extinction_tolerance
import ResilienceCore: jacobian, SignPattern, quasi_sign_stability, press_sign_determinacy

export reactivity, resilience, amplification_envelope, max_amplification, pulse_covariance, median_reactivity,
       median_displacement, median_max_amplification, median_resilience, stochastic_invariability,
       sensitivity_matrix, press_sensitivity, total_sensitivity
export PerturbationOptions, global_mortality_tolerance, local_mortality_tolerances, press_response,
       removal_responses, extinction_tolerance, jacobian
export RAW_COLUMNS, METRICS, community_record, published_pulse_covariance, amplification_metrics, prepare_metrics,
       read_metrics
export MatrixMetadata, InteractionMatrix, matrix, matrix_metadata, taxa, nspecies, AlphaMatrix, glv_rhs!, equilibrium,
       isfeasible, Jacobian, CommunityMatrix, per_capita, FluxMatrix, SignedDigraph, ProbabilisticNetwork,
       jacobian_metrics, subcommunity, invasion_rates, steady_state, glv_trajectory
export feasibility_domain, structural_vector, feasibility_angle, quasi_sign_stability, press_sign_determinacy,
       SignPattern
export smap_jacobians, invariability
export read_kinlock, alpha_from_rii, read_adler, read_stein
export SizeCorrelations, size_correlations, class_average, SIZE_CLASSES, metric_network, modularity,
       modular_groups, cophenetic, metric_dendrogram, explained_variance, ellipsoid_volume, correlation_types

include("local_metrics.jl")
include("dataset.jl")
include("analysis.jl")
include("matrices.jl")
include("structural.jl")
include("timeseries.jl")
include("readers.jl")

end
