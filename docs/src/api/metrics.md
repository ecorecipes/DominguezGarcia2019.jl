# Datasets and conventions

The stability metrics themselves are ResilienceCore's, re-exported here: `resilience`,
`reactivity`, `amplification_envelope`, `max_amplification`, `pulse_covariance`, the median
metrics, `stochastic_invariability`, the press sensitivities, and the simulated experiments
(`PerturbationOptions`, `global_mortality_tolerance`, `local_mortality_tolerances`,
`press_response`, `removal_responses`, `extinction_tolerance`). Their documentation is in
[ResilienceCore.jl](https://ecorecipes.github.io/ResilienceCore.jl/dev/). The experiments run
on the assembled communities of
[FoodWebDynamics.jl](https://ecorecipes.github.io/FoodWebDynamics.jl/dev/), which implement
ResilienceCore's community interface.

This page documents what is specific to the paper: the Jacobian of an assembled community,
the conventions the published Dataset S1 actually used, and the dataset layout.

## The community's Jacobian

```@docs
community_jacobian
```

## The published conventions

```@docs
published_max_amplification
published_pulse_covariance
amplification_metrics
```

## Datasets

```@docs
RAW_COLUMNS
METRICS
community_record
prepare_metrics
read_metrics
```
