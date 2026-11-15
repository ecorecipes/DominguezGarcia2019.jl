# EcologicalResilience.jl

Resilience and stability metrics for ecological models and data. The first target is an
end-to-end reproduction of Domínguez-García, Dakos & Kéfi (2019), "Unveiling dimensions of
stability in complex ecological networks", *PNAS* 116:25714
([doi:10.1073/pnas.1904470116](https://doi.org/10.1073/pnas.1904470116)).

- **Food webs**: the niche model of Williams & Martinez (2000) and an allometric bioenergetic
  model, with community assembly as in SI §8–9. The authors' simulator is not public; the
  details the paper leaves open were settled with their published dataset.
- **Metrics**: the paper's 27 stability metrics. Some come from the Jacobian (reactivity,
  amplification, resilience, their median versions, stochastic invariability, sensitivity);
  the others come from simulated presses and extinctions (tolerance, resistance, sensitivity,
  cascading extinctions, robustness). The Jacobian-based functions take any Jacobian, not only
  the food-web model's.
- **Analysis**: per-size Spearman correlations, size-class averages, networks of metrics and
  their modular groups, hierarchical clustering, explained variance and covariance ellipsoids.
  The results are checked against the authors' Python code.
- **Interaction matrices and data**: typed interaction matrices, Lotka–Volterra assembly,
  structural stability and sign-pattern analyses. Readers cover empirical datasets (Kinlock
  plant interactions, Adler and Stein communities) and observed food webs (GATEWAy, aggregated
  into trophic species).

```@docs
EcologicalResilience
```
