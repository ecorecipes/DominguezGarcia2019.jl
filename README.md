# EcologicalResilience.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://ecorecipes.github.io/EcologicalResilience.jl/dev/)

Resilience and stability metrics for ecological models and data, starting with an end-to-end
reproduction of Domínguez-García, Dakos & Kéfi (2019), *PNAS* 116:25714, "Unveiling dimensions
of stability in complex ecological networks".

It consumes models from other ecorecipes packages rather than defining them:

- **Metrics on operators** come from
  [ResilienceCore.jl](https://github.com/ecorecipes/ResilienceCore.jl) and are re-exported
  (`resilience`, `reactivity`, `max_amplification`, `stochastic_invariability`,
  `median_resilience`, the simulated presses and extinctions, `applicable_metrics`, …).
- **Simulation** lives in [FoodWebDynamics.jl](https://github.com/ecorecipes/FoodWebDynamics.jl)
  (`niche_model`, `FoodWebModel`, `assemble_community`, ATN variants). Loading it activates the
  `EcologicalResilienceFoodWebDynamicsExt` extension, which adds `community_record` (all 27
  metrics, with the raw columns of Dataset S1) and `Jacobian(c)` for assembled communities. The
  authors' simulator is not public; the defaults reproduce summaries of their Dataset S1, which
  also settles several points where the SI is silent or wrong.
- **Analysis** (`size_correlations`, `class_average`, `modular_groups`, `metric_dendrogram`,
  `explained_variance`, `ellipsoid_volume`, `correlation_types`): the paper's analysis, checked
  against the authors' Python code on their Dataset S2. It needs only a metrics table.

```julia
using EcologicalResilience, FoodWebDynamics, Random

rng = Xoshiro(1)
c = nothing
while c === nothing
    global c, _ = assemble_community(rng, 30)
end
J = Jacobian(c)                         # the Jacobian at the assembled equilibrium
resilience(J), reactivity(J), max_amplification(J)
rec = community_record(rng, c)          # all raw columns of Dataset S1
```

Tutorials (Quarto, in `vignettes/`): reproducing the analysis on the authors' data; one food
web and its 27 metrics; the simulated ensemble compared with the published one; the
early-response group as an artefact of the published pulse weights.

`data/dominguez-garcia2019/` holds the authors' Dataset S2 and one of their group files, from
<https://github.com/domgarvir/Stability-metrics> (no licence stated there).
`data/simulated/ensemble_v4.csv` is 10,000 communities simulated with
`scripts/simulate_dataset.jl` (seed 20191217, initial richness 5–110), used by the third
tutorial; each row can be regenerated from its `network` seed. Its raw columns follow the
published convention of Dataset S1 (`community_record`); `*_si` columns hold the SI's
versions of `Amax`, `tmax` and the median metrics. Design notes and
the verification record are in the workspace wiki (`../../wiki`).

## Status

Not registered; the ecorecipes dependencies are fetched by URL through `[sources]`
(Julia ≥ 1.11). Tests: `julia --project -e 'using Pkg; Pkg.test()'`.
