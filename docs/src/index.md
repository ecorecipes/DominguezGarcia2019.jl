# DominguezGarcia2019.jl

A reproduction of Domínguez-García, Dakos & Kéfi (2019), "Unveiling dimensions of stability
in complex ecological networks", *PNAS* 116:25714
([doi:10.1073/pnas.1904470116](https://doi.org/10.1073/pnas.1904470116)).

The paper computes 27 stability metrics on simulated niche-model food webs, correlates them
within each species richness, and groups them into three "dimensions" of stability. This
package reproduces it end to end, from the simulator, which is not public and was rebuilt from
the SI, to the published groups, and then asks how far the conclusions carry.

It is a case study. It holds only what is specific to the paper, and takes everything else from
the ecorecipes packages:

| Package | Provides |
|---|---|
| [ResilienceCore.jl](https://ecorecipes.github.io/ResilienceCore.jl/dev/) | the metrics (resilience, reactivity, amplification, median metrics, invariability, sensitivities) and the press and extinction experiments; re-exported here |
| [FoodWebDynamics.jl](https://ecorecipes.github.io/FoodWebDynamics.jl/dev/) | the niche model, the allometric bioenergetic model, community assembly and the GATEWAy reader |
| [CommunityMatrices.jl](https://ecorecipes.github.io/CommunityMatrices.jl/dev/) | interaction matrices, structural stability and sign patterns, which began in this package |
| [DynamicalResilience.jl](https://ecorecipes.github.io/DynamicalResilience.jl/dev/) | resilience of SciML models: steady states, thresholds, recovery and basin stability |

What is here:

- **Conventions:** [`community_record`](@ref) computes the 39 raw columns of Dataset S1 for an
  assembled community. Two of them follow what the published data did rather than what the
  SI says: the ∞-norm for `Amax` and `tmax`, and misaligned pulse weights for the median
  metrics. Both can be switched to the SI's definitions.
- **Dataset transforms:** [`prepare_metrics`](@ref) and [`read_metrics`](@ref) turn raw
  columns into the paper's 27 "larger is more stable" metrics, as the authors' code does.
- **Analysis:** per-size Spearman correlations, size-class averages, networks of metrics and
  their modular groups, hierarchical clustering, explained variance and covariance ellipsoids.
  They are checked against the authors' Python code.
- **Data:** the authors' Dataset S2, a 10,000-community ensemble simulated with the
  reconstructed model, and ensembles under alternative model settings.

The tutorials reproduce the analysis on the authors' data, follow one food web through all
27 metrics, compare a simulated ensemble with the published one, show that the "early response
to pulse" dimension is an indexing artefact, run the pipeline on observed food webs, and test
the conclusions against changes to the model's settings.

```@docs
DominguezGarcia2019
```
