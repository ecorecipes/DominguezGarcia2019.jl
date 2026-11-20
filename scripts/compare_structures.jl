# Compare the correlation structure of the 27 stability metrics between datasets: the published
# Dataset S2, a niche-model ensemble simulated with this package, and observed food-web
# topologies (scripts/observed_webs.jl).
#
#   julia --project scripts/compare_structures.jl PUBLISHED_S2.csv SIMULATED.csv OBSERVED.csv
#
# Observed webs are few at any one richness, so every dataset is binned identically by final
# richness (5–15, 16–30, 31–60, 61–120); Spearman correlations are computed within bins and
# averaged over bins. Rows with a status other than "ok", or any non-finite metric, are dropped.
#
# Differences between datasets are judged against a noise ceiling: the same comparison between
# two random halves of one dataset (20 splits). Observed communities are split by web, so the
# replicates of a web stay together.
#
# Metrics that are constant in any dataset (their correlations are undefined) are dropped from
# all datasets, so every comparison uses the same set of metrics.
#
# The simulated ensemble follows the published convention for Amax, tmax and the median metrics
# (`community_record`). If it also has the SI versions (`*_si`, from simulate_dataset.jl), a
# fourth dataset `niche_si` uses them, and it is the one compared with observed webs, which are
# run with the SI convention.
using DominguezGarcia2019
using CSV, DataFrames, Statistics, Printf
using Random: Xoshiro, shuffle

const BINS = [(5, 15), (16, 30), (31, 60), (61, 120)]

const SI_COLUMNS = [:Amax, :tmax, :MRinf, :MR0, :MAmax, :Mtmax]

# si = true: replace the published-convention columns by the SI versions (`*_si`), which observed
# webs use (they have no niche order); see `amplification_metrics`
function load(path; si = false)
    raw = CSV.read(path, DataFrame)
    hasproperty(raw, :status) && (raw = raw[coalesce.(raw.status .== "ok", false), :])
    si && for c in SI_COLUMNS
        raw[!, c] = raw[!, Symbol(c, "_si")]
    end
    # observed webs are pressed as the SI says (10%); so is niche_si where the ensemble has it
    si && hasproperty(raw, :RM_p_si) && (raw.RM_p = raw.RM_p_si; raw.L1DBRM_p = raw.L1DBRM_p_si)
    unit = hasproperty(raw, :web) ? string.(raw.web) : string.(1:nrow(raw))
    m = hasproperty(raw, :L1DBRMt_s) ? prepare_metrics(raw) : read_metrics(path)
    m.unit = unit                         # the sampling unit for split halves
    m = m[[all(isfinite, m[i, METRICS]) for i in 1:nrow(m)], :]
    b = [something(findfirst(r -> r[1] ≤ n ≤ r[2], BINS), 0) for n in m.N]
    m = m[b .> 0, :]
    m.N = b[b .> 0]                       # size_correlations groups by N: use the bin index
    return m
end

# mean over bins with at least `min_n` communities, ignoring correlations undefined because a
# metric is constant within a bin
function binned_correlations(m; min_n = 20)
    sc = size_correlations(m; labels = LABELS, sizes = 1:length(BINS), max_sample = typemax(Int))
    use = findall(≥(min_n), sc.n)
    k = length(LABELS)
    R = [let v = filter(!isnan, sc.rho[i, j, use]); isempty(v) ? 0.0 : mean(v) end for i in 1:k, j in 1:k]
    return R, sc, use
end

function structure(name, m)
    R, sc, use = binned_correlations(m)
    for k in use, (j, l) in enumerate(LABELS)
        all(isnan, (sc.rho[j, i, k] for i in eachindex(LABELS) if i != j)) &&
            println("    $name: $l is constant in bin $(BINS[sc.sizes[k]])")
    end
    g, Q = modular_groups(metric_network(R); runs = 50, rng = Xoshiro(1))
    @printf("%-10s communities per bin %s (bins used: %d); Louvain Q %.3f, %d groups\n", name, sc.n, length(use), Q,
            maximum(g))
    for k in 1:maximum(g)
        println("    ", join(LABELS[g .== k], " "))
    end
    return (R = R, groups = g)
end

upper(R) = [R[i, j] for j in axes(R, 2) for i in axes(R, 1) if i < j]
rand_index(a, b) = mean((a[i] == a[j]) == (b[i] == b[j]) for i in eachindex(a) for j in eachindex(a) if i < j)

function split_half(m; splits = 20)
    units = unique(m.unit)
    r, d = Float64[], Float64[]
    for s in 1:splits
        half = Set(shuffle(Xoshiro(s), units)[1:(length(units) ÷ 2)])
        inhalf = in.(m.unit, Ref(half))
        x, y = upper(binned_correlations(m[inhalf, :])[1]), upper(binned_correlations(m[.!inhalf, :])[1])
        push!(r, cor(x, y))
        push!(d, mean(abs.(x .- y)))
    end
    return r, d
end

names_ = ["published", "niche", "observed"]
data = Dict(name => load(path) for (name, path) in zip(names_, ARGS))
# niche webs against observed webs: both with the SI's definitions, when the ensemble has them
nicheobs = "niche"
if hasproperty(CSV.read(ARGS[2], DataFrame; limit = 1), :Amax_si)
    data["niche_si"] = load(ARGS[2]; si = true)
    push!(names_, "niche_si")
    nicheobs = "niche_si"
end
constant = [l for l in METRICS if any(length(unique(data[n][!, l])) == 1 for n in names_)]
const LABELS = setdiff(METRICS, constant)
isempty(constant) || println("Constant in some dataset, dropped from all: ", join(constant, ", "), "\n")
res = Dict(name => structure(name, data[name]) for name in names_)

println("\nBetween datasets")
for (a, b) in (("published", "niche"), ("published", "observed"), (nicheobs, "observed"))
    x, y = upper(res[a].R), upper(res[b].R)
    @printf("  %-9s vs %-9s: correlation of pair correlations %.2f, mean |Δρ| %.3f, Rand index of groups %.2f\n", a, b,
            cor(x, y), mean(abs.(x .- y)), rand_index(res[a].groups, res[b].groups))
end

println("\nNoise ceiling: two random halves of one dataset (20 splits)")
for name in names_
    r, d = split_half(data[name])
    @printf("  %-9s: correlation of pair correlations %.2f (range %.2f–%.2f), mean |Δρ| %.3f\n", name, mean(r),
            extrema(r)..., mean(d))
end

println("\nMetric pairs that differ most between niche-model and observed webs")
pairs = [(i, j) for j in eachindex(LABELS) for i in eachindex(LABELS) if i < j]
Δ = [res["observed"].R[i, j] - res[nicheobs].R[i, j] for (i, j) in pairs]
for k in sortperm(abs.(Δ); rev = true)[1:10]
    i, j = pairs[k]
    @printf("  %-6s %-6s published %5.2f  %s %5.2f  observed %5.2f\n", LABELS[i], LABELS[j],
            res["published"].R[i, j], nicheobs, res[nicheobs].R[i, j], res["observed"].R[i, j])
end
