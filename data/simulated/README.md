# Simulated ensemble

`ensemble_v4.csv` holds 10,000 niche-model communities simulated with this package:

```bash
julia -t 10 --gcthreads=1 --project scripts/simulate_dataset.jl data/simulated/ensemble_v4.csv 10000
```

The defaults are seed 20191217 and initial richness uniform on 5–110. Community `Sim<seed + k>`
uses the random stream `Xoshiro(seed + k)`, so every row can be regenerated.

Columns:
- `network`;
- the 39 raw columns of Dataset S1 (`RAW_COLUMNS`), following the published convention
  (`community_record(...; convention = :published)`);
- `attempts`, the number of assembly draws until a viable community;
- `seconds`, the time spent on its metrics;
- `Amax_si`, `tmax_si`, `MRinf_si`, `MR0_si`, `MAmax_si` and `Mtmax_si`, the same six metrics
  as the SI defines them (`amplification_metrics(c; convention = :si)`);
- `RM_p_si` and `L1DBRM_p_si`, the global press as the SI defines it (every mortality raised by
  10%, `global_press_response`). `RM_p` and `L1DBRM_p` hold the published press, an absolute
  rise of 0.055 in every mortality (`published_press_factors`).

The ensemble was first simulated with the 10% press in `RM_p` and `L1DBRM_p`. Those two columns
were recomputed with the published press, and the 10% values moved to the `_si` columns, by

```bash
julia -t 10 --gcthreads=1 --project scripts/recompute_global_press.jl data/simulated/ensemble_v4.csv PRESS.csv
```

which regenerates each community from its seed and checks it against the stored `Rinf` (all
10,000 matched). `simulate_dataset.jl` now writes both presses directly.

It was generated with the package defaults: no plant metabolism, metabolic exponent −0.25,
h₀ = 0.4, the survival rule `N ≥ floor(0.9 iniN)`, and diet weights kept after extinctions.

## Robustness variants

`robustness/` holds six ensembles that each change one model setting, on the same community
seeds 1–1,200 (`Sim20191218` to `Sim20192417`) as the first 1,200 rows of `ensemble_v4.csv`.
The vignette *Robustness of the reproduction to model settings* analyses them. Their columns
are those of `ensemble_v4.csv`. They were simulated before the published press was adopted;
their `RM_p` and `L1DBRM_p` were then recomputed with it, and the 10% values moved to
`RM_p_si` and `L1DBRM_p_si`. The recompute takes the variant's setting as the simulation did:

```bash
julia -t 10 --gcthreads=1 --project scripts/recompute_global_press.jl \
    data/simulated/robustness/h0_0.2.csv OUT.csv 5 110 h0=0.2
```

Every community re-assembled exactly, with the Rinf stored in its row.

| File | Setting (`KEY=VALUE`) | Communities |
|---|---|---|
| `h0_0.2.csv` | `h0=0.2` | 1,200 |
| `h0_0.7.csv` | `h0=0.7` | 95 |
| `q_0.6.csv` | `q=0.6` | 1,200 |
| `q_0.1.csv` | `q=0.1` | 6 |
| `plant_metabolism.csv` | `plant_metabolism=true` | 1,194 |
| `metabolic_exponent_+0.25.csv` | `metabolic_exponent=0.25` | 9 |

Each was simulated with

```bash
julia -t 10 --gcthreads=1 --project scripts/simulate_dataset.jl OUT.csv 1200 5 110 20191217 KEY=VALUE
```

The `h0=0.7`, `q=0.1` and `metabolic_exponent=0.25` runs rarely reach a viable community and
were stopped after a 90-minute budget, so they hold only the communities reached by then.
`plant_metabolism=true` holds 1,194 of the 1,200 seeds. To print the summary table:

```bash
julia --project scripts/robustness.jl data/dominguez-garcia2019/Dataset_S2.csv 1200 \
    baseline=data/simulated/ensemble_v4.csv h0_0.2=data/simulated/robustness/h0_0.2.csv ...
```
