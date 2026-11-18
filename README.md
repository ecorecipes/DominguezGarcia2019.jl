# DominguezGarcia2019.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://ecorecipes.github.io/DominguezGarcia2019.jl/dev/)

An end-to-end reproduction of Domínguez-García, Dakos & Kéfi (2019), *PNAS* 116:25714,
"Unveiling dimensions of stability in complex ecological networks": 27 stability metrics on
simulated food webs, their correlations within each species richness, and the three
"dimensions" of stability they form.

This is a case study; it holds only what is specific to the paper.

- **Metrics** come from [ResilienceCore.jl](https://github.com/ecorecipes/ResilienceCore.jl)
  and are re-exported: `resilience`, `reactivity`, `max_amplification`, the median metrics,
  `stochastic_invariability`, the sensitivities, and the simulated presses and extinctions.
- **Simulation** comes from [FoodWebDynamics.jl](https://github.com/ecorecipes/FoodWebDynamics.jl):
  `niche_model`, `FoodWebModel`, `assemble_community`, `read_gateway`. The authors' simulator
  is not public; the defaults reproduce summaries of their Dataset S1, which also settles
  several points where the SI is silent or wrong.
- **Here:**
  - `community_record` gives all 39 raw columns of Dataset S1 for an assembled community, with
    the conventions the published data used (`published_max_amplification`,
    `published_pulse_covariance`) or the SI's (`convention = :si`).
  - `prepare_metrics` and `read_metrics` turn raw columns into the paper's 27 metrics.
  - `size_correlations`, `class_average`, `modular_groups`, `metric_dendrogram`,
    `explained_variance`, `ellipsoid_volume` and `correlation_types` port the paper's
    analysis, checked against the authors' Python code on their Dataset S2.

Interaction matrices, structural stability and sign patterns, which began in this package
(as EcologicalResilience.jl), are now in
[CommunityMatrices.jl](https://github.com/ecorecipes/CommunityMatrices.jl).

```julia
using DominguezGarcia2019, FoodWebDynamics, Random

rng = Xoshiro(1)
c = nothing
while c === nothing
    global c, _ = assemble_community(rng, 30)
end
op = community_jacobian(c)              # ResilienceCore's ContinuousJacobian at the equilibrium
resilience(op), reactivity(op), max_amplification(op)
rec = community_record(rng, c)          # all raw columns of Dataset S1
```

Tutorials (Quarto, in `vignettes/`):
1. reproducing the analysis on the authors' data;
2. one food web and its 27 metrics;
3. the simulated ensemble compared with the published one;
4. the early-response group as an artefact of the published pulse weights;
5. observed food webs (GATEWAy, downloaded at render time to
   `~/.cache/DominguezGarcia2019/gateway`);
6. robustness of the conclusions to the model's settings.

`data/dominguez-garcia2019/` holds the authors' Dataset S2 and one of their group files, from
<https://github.com/domgarvir/Stability-metrics> (no licence stated there).
`data/simulated/ensemble_v4.csv` is 10,000 communities simulated with
`scripts/simulate_dataset.jl` (seed 20191217, initial richness 5–110); each row can be
regenerated from its `network` seed. Its raw columns follow the published convention of
Dataset S1 (`community_record`); `*_si` columns hold the SI's versions of `Amax`, `tmax` and
the median metrics. `data/simulated/robustness/` holds the ensembles under alternative model
settings (see `data/simulated/README.md`).

## Status

Not registered; the ecorecipes dependencies are fetched by URL through `[sources]`
(Julia ≥ 1.11). Tests: `julia --project -e 'using Pkg; Pkg.test()'`.
