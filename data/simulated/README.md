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
  as the SI defines them (`amplification_metrics(c; convention = :si)`).

It was generated with the package defaults: no plant metabolism, metabolic exponent −0.25,
h₀ = 0.4, the survival rule `N ≥ floor(0.9 iniN)`, and diet weights kept after extinctions.
