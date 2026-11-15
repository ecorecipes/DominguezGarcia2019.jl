# Domínguez-García, Dakos & Kéfi (2019): Dataset S2 and group file

Copied unchanged from the authors' repository
[domgarvir/Stability-metrics](https://github.com/domgarvir/Stability-metrics), commit
`b1cc405`, directory `DATA/`. The data accompany:

> Domínguez-García V, Dakos V, Kéfi S (2019) Unveiling dimensions of stability in complex
> ecological networks. *PNAS* 116:25714.
> [doi:10.1073/pnas.1904470116](https://doi.org/10.1073/pnas.1904470116)

| file | content |
|---|---|
| `Dataset_S2.csv` | the paper's analysed sample: 100 simulated communities per final richness (N 5–99), raw metric columns as in the authors' Dataset S1 (`network`, 28 metric columns, `N`, `iniN`, `Nb`, `Napex`, `Ncan`, `Btot`, `B_IPR`, `a_i`, `a_j`, `h_i`, `h_j`) |
| `groups_rMED.graphml` | the authors' Gephi partition of the metric network into stability groups (their three size-class files are byte-identical) |

The raw columns are not yet sign-adjusted. `read_metrics` applies the authors'
`prepare_dataframe` transformations; the publisher's version of Dataset S2 (PNAS SI) holds the
transformed values.

**Licence.** The source repository states no licence. The files are redistributed here only
so that the tests and vignettes can reproduce the paper's analysis; the authors keep all
rights. Please cite the paper above when using them.
