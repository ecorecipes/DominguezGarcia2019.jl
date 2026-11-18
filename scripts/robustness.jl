# How robust are the paper's stability dimensions to modelling choices? Each simulated
# ensemble (scripts/simulate_dataset.jl, possibly with KEY=VALUE model settings) is analysed
# under both median-metric conventions and compared with the published Dataset S2.
#
#   julia --project scripts/robustness.jl PUBLISHED_S2.csv NSEEDS NAME=ENSEMBLE.csv ...
#
# Only the communities with seeds 1..NSEEDS (network Sim20191218 ...) are used from each
# ensemble, so every variant is compared on the same draws. Correlations are computed within
# richness bins (5–15, 16–30, 31–60, 61–99) and averaged, as in compare_structures.jl.
using DominguezGarcia2019
using CSV, DataFrames, Statistics, Printf
using Random: Xoshiro

const BINS = [(5, 15), (16, 30), (31, 60), (61, 99)]
const SEED0 = 20191217
const SIX = (:Amax, :tmax, :MRinf, :MR0, :MAmax, :Mtmax)

function binned(m)
    b = [something(findfirst(r -> r[1] ≤ n ≤ r[2], BINS), 0) for n in m.N]
    m = m[b .> 0, :]
    m.N = b[b .> 0]
    sc = size_correlations(m; sizes = 1:length(BINS), max_sample = typemax(Int))
    k = length(sc.labels)
    R = [let v = filter(!isnan, sc.rho[i, j, :]); isempty(v) ? NaN : mean(v) end for i in 1:k, j in 1:k]
    return R, sc.labels
end

function structure(m)
    R, labels = binned(m)
    g, Q = modular_groups(metric_network(R); runs = 50, rng = Xoshiro(1))
    return (R = R, labels = labels, groups = g, Q = Q)
end

upper(R) = [R[i, j] for j in axes(R, 2) for i in axes(R, 1) if i < j]
rand_index(a, b) = mean((a[i] == a[j]) == (b[i] == b[j]) for i in eachindex(a) for j in eachindex(a) if i < j)
function early_group(s)                       # MR0, MAmax and Mtmax in one group, and nothing constant
    i = [findfirst(==(l), s.labels) for l in ("MR0", "MAmax", "Mtmax")]
    return length(unique(s.groups[i])) == 1 ? "yes" : "no"
end

pubpath, nseeds = ARGS[1], parse(Int, ARGS[2])
ref = structure(read_metrics(pubpath))
println("published: ", maximum(ref.groups), " groups, Q ", round(ref.Q; digits = 3))
for k in 1:maximum(ref.groups)
    println("    ", join(ref.labels[ref.groups .== k], " "))
end
@printf("\n%-14s %-10s %5s %6s %6s %7s %6s %6s %8s %7s %6s %6s\n", "variant", "weights", "n", "attempts", "N/iniN",
        "Nb/N", "MAmax>1", "r(ρ)", "mean|Δρ|", "Rand", "groups", "early")
for arg in ARGS[3:end]
    name, path = split(arg, "="; limit = 2)
    raw = CSV.read(path, DataFrame)
    raw = raw[in.(raw.network, Ref(Set("Sim$(SEED0 + k)" for k in 1:nseeds))) .& (5 .<= raw.N .<= 99), :]
    for weights in ("published", "SI")
        d = copy(raw)
        weights == "SI" && foreach(c -> (d[!, c] = d[!, Symbol(c, "_si")]), SIX)
        s = structure(prepare_metrics(d))
        x, y = upper(ref.R), upper(s.R)
        ok = isfinite.(x) .& isfinite.(y)
        @printf("%-14s %-10s %5d %8.1f %6.3f %7.3f %6.1f%% %6.2f %8.3f %7.2f %6d %6s\n", name, weights, nrow(d),
                mean(d.attempts), median(d.N ./ d.iniN), median(d.Nb ./ d.N), 100mean(d.MAmax .> 1 + 1e-12),
                cor(x[ok], y[ok]), mean(abs.(x[ok] .- y[ok])), rand_index(ref.groups, s.groups), maximum(s.groups),
                early_group(s))
    end
end
