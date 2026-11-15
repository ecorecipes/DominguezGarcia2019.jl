# Interaction matrices

Typed matrices make explicit which object a number belongs to, and conversions between them
ask for what they need.

```@docs
MatrixMetadata
InteractionMatrix
matrix
matrix_metadata
taxa
nspecies
AlphaMatrix
glv_rhs!
equilibrium
isfeasible
Jacobian
CommunityMatrix
per_capita
FluxMatrix
SignedDigraph
ProbabilisticNetwork
jacobian_metrics
```

## Lotka–Volterra assembly and trajectories

```@docs
subcommunity
invasion_rates
steady_state
glv_trajectory
```

## Structural stability and sign patterns

```@docs
feasibility_domain
structural_vector
feasibility_angle
quasi_sign_stability
press_sign_determinacy
```

## From time series

```@docs
smap_jacobians
invariability
```

## Empirical data

Each reader records its conversions in the matrices' metadata. The data and their provenance
are described in `lab/data/*/README.md`.

```@docs
read_kinlock
alpha_from_rii
read_adler
read_stein
ObservedWeb
read_gateway
trophic_species
relative_masses
```

## Symbolic models

```@docs
reaction_network
```
