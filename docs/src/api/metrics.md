# Stability metrics

## From the Jacobian

```@docs
reactivity
resilience
amplification_envelope
max_amplification
pulse_covariance
median_reactivity
median_displacement
median_max_amplification
median_resilience
stochastic_invariability
sensitivity_matrix
press_sensitivity
total_sensitivity
```

## From simulated press perturbations

The experiments `PerturbationOptions`, `global_mortality_tolerance`,
`local_mortality_tolerances`, `press_response`, `removal_responses` and
`extinction_tolerance` are defined in
[ResilienceCore.jl](https://ecorecipes.github.io/ResilienceCore.jl/dev/) and re-exported here.
They run on any community that implements ResilienceCore's interface (`state`,
`present_species`, `with_mortality`, `without_species`, `jacobian`), such as the assembled
communities of FoodWebDynamics.jl.

## Datasets

```@docs
RAW_COLUMNS
METRICS
community_record
amplification_metrics
published_pulse_covariance
prepare_metrics
read_metrics
```
