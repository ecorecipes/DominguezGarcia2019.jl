# EcologicalResilience.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://ecorecipes.github.io/EcologicalResilience.jl/dev/)

Resilience and stability metrics for ecological models and data, starting with an end-to-end
reproduction of Domínguez-García, Dakos & Kéfi (2019), *PNAS* 116:25714, "Unveiling dimensions
of stability in complex ecological networks".

- **Simulation** (`niche_model`, `FoodWebModel`, `assemble_community`): niche-model food webs
  with allometric bioenergetic dynamics, assembled to a locally stable equilibrium as in SI §8–9.
  The authors' simulator is not public. The defaults reproduce summaries of their Dataset S1,
  which also settles several points where the SI is silent or wrong.
- **Metrics** (`community_record`): all 27 metrics, with the raw columns of Dataset S1. The
  Jacobian metrics (`reactivity`, `resilience`, `max_amplification`, `stochastic_invariability`,
  `median_resilience`, …) work on any Jacobian.
- **Analysis** (`size_correlations`, `class_average`, `modular_groups`, `metric_dendrogram`,
  `explained_variance`, `ellipsoid_volume`, `correlation_types`): the paper's analysis, checked
  against the authors' Python code on their Dataset S2.

```julia
using EcologicalResilience, Random

rng = Xoshiro(1)
c = nothing
while c === nothing
    global c, _ = assemble_community(rng, 30)
end
J = jacobian(c.B, c.model)
resilience(J), reactivity(J), max_amplification(J)
rec = community_record(rng, c)          # all raw columns of Dataset S1
```

Tutorials (Quarto, in `vignettes/`): reproducing the analysis on the authors' data; one food
web and its 27 metrics; the simulated ensemble compared with the published one.

`data/dominguez-garcia2019/` holds the authors' Dataset S2 and one of their group files, from
<https://github.com/domgarvir/Stability-metrics> (no licence stated there).
`data/simulated/ensemble_v4.csv` is 3,000 communities simulated with
`scripts/simulate_dataset.jl` (seed 20191217, initial richness 5–110), used by the third
tutorial; each row can be regenerated from its `network` seed. Its raw columns follow the
published convention of Dataset S1 (`community_record`); `*_si` columns hold the SI's
versions of `Amax`, `tmax` and the median metrics. Design notes and
the verification record are in the workspace wiki (`../../wiki`).

## Status

Not registered. Tests: `julia --project -e 'using Pkg; Pkg.test()'`.
