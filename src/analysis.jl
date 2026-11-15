# The statistical analysis of Domínguez-García et al. (2019): Spearman correlations of the
# stability metrics at each community size, their averages over size classes, the network of
# metrics and its modular groups, hierarchical clustering, explained variance and the volume of
# covariance ellipsoids. Behaviour follows the authors' Functions.py; where the code and the
# paper differ, both are available (see lab/wiki/analysis/).

"""
    SizeCorrelations

Per-size correlation structure of a metrics table: for each richness in `sizes`, the Spearman
correlation matrix `rho`, the covariance of ranks `rankcov` and the covariance of
mean/range-normalised values `valuecov` of the metrics in `labels` (arrays are
`metrics × metrics × sizes`), and the number of communities `n` used at each size.
"""
struct SizeCorrelations
    labels::Vector{String}
    sizes::Vector{Int}
    n::Vector{Int}
    rho::Array{Float64,3}
    rankcov::Array{Float64,3}
    valuecov::Array{Float64,3}
end

rangenormalize(x) = (x .- mean(x)) ./ (maximum(x) - minimum(x))

"""
    size_correlations(df; labels = METRICS, sizes = 5:99, max_sample = 100, rng = Random.default_rng())

Correlations among the metrics `labels` of a table such as [`prepare_metrics`](@ref) returns,
separately for the communities of each richness `N` in `sizes`, as in the authors'
`get_full_correlation_bysize` (window width 0). At most `max_sample` communities per size are
used, drawn with `rng` when there are more. The authors' loop covers N = 5–99.
"""
function size_correlations(df::DataFrame; labels::AbstractVector{<:AbstractString} = METRICS,
                           sizes = 5:99, max_sample::Integer = 100, rng::AbstractRNG = Random.default_rng())
    k = length(labels)
    used = Int[]
    rho = Array{Float64,3}(undef, k, k, 0)
    rc = similar(rho)
    vc = similar(rho)
    ns = Int[]
    for s in sizes
        rows = findall(==(s), df.N)
        length(rows) < 3 && continue
        length(rows) > max_sample && (rows = StatsBase.sample(rng, rows, max_sample; replace = false))
        X = Matrix{Float64}(df[rows, labels])
        R = reduce(hcat, (StatsBase.tiedrank(c) for c in eachcol(X)))
        V = reduce(hcat, (rangenormalize(c) for c in eachcol(X)))
        rho = cat(rho, cor(R); dims = 3)
        rc = cat(rc, cov(R); dims = 3)
        vc = cat(vc, cov(V); dims = 3)
        push!(used, s)
        push!(ns, length(rows))
    end
    return SizeCorrelations(collect(String, labels), used, ns, rho, rc, vc)
end

sizeindex(sc::SizeCorrelations, lo::Integer, hi::Integer) = findall(s -> lo ≤ s ≤ hi, sc.sizes)

"""
    class_average(sc::SizeCorrelations, lo, hi; field = :rho)

The element-wise mean and standard deviation over the sizes `lo:hi` (inclusive) of one of the
per-size matrices (`:rho`, `:rankcov` or `:valuecov`), as the authors' `get_av_corr_from_sdf`
and `get_av_cov_from_sdf` compute for the size classes small (5–15), medium (45–55) and large
(85–95). Like pandas' `mean()` and `std()` there, undefined values (`NaN`, for a metric that
is constant at some size) are skipped; an entry undefined at every size stays `NaN`.
"""
function class_average(sc::SizeCorrelations, lo::Integer, hi::Integer; field::Symbol = :rho)
    A = getfield(sc, field)[:, :, sizeindex(sc, lo, hi)]
    skipnan(f, v) = (x = filter(!isnan, v); isempty(x) ? NaN : length(x) == 1 && f === std ? NaN : f(x))
    return ([skipnan(mean, A[i, j, :]) for i in axes(A, 1), j in axes(A, 2)],
            [skipnan(std, A[i, j, :]) for i in axes(A, 1), j in axes(A, 2)])
end

"""
    SIZE_CLASSES

The three size classes of the paper, as `name => (lo, hi)` in the authors' processing order.
"""
const SIZE_CLASSES = ["small" => (5, 15), "large" => (85, 95), "medium" => (45, 55)]

"""
    metric_network(R)

The weighted, unsigned network of metrics from a matrix of mean correlations `R`: edge weights
`|Rᵢⱼ|` (the absolute value of the averaged correlation, as in the authors' code), no
self-loops. An undefined correlation (`NaN`) gives no edge, as the authors' threshold test
`abs(mean) > threshold` is false for `NaN`.
"""
function metric_network(R::AbstractMatrix)
    W = [isnan(r) ? 0.0 : abs(r) for r in R]
    W[diagind(W)] .= 0
    return W
end

"""
    modularity(W, groups; resolution = 1.0)

Newman's modularity of the partition `groups` (one label per node) of the weighted undirected
network `W`.
"""
function modularity(W::AbstractMatrix, groups::AbstractVector; resolution::Real = 1.0)
    k = vec(sum(W; dims = 2))
    twom = sum(k)
    Q = 0.0
    for i in axes(W, 1), j in axes(W, 2)
        groups[i] == groups[j] && (Q += W[i, j] - resolution * k[i] * k[j] / twom)
    end
    return Q / twom
end

# One run of the Louvain method (Blondel et al. 2008): local moves in random order, then
# aggregation, repeated until nothing changes. Returns a label per original node.
function louvain_once(rng::AbstractRNG, W::AbstractMatrix, resolution::Real)
    n = size(W, 1)
    labels = collect(1:n)
    A = Matrix{Float64}(W)
    members = [[i] for i in 1:n]
    while true
        m = size(A, 1)
        k = vec(sum(A; dims = 2))
        twom = sum(k)
        comm = collect(1:m)
        tot = copy(k)
        improved = false
        moved = true
        while moved
            moved = false
            for i in Random.shuffle(rng, 1:m)
                ci = comm[i]
                links = Dict{Int,Float64}()
                for j in 1:m
                    j == i && continue
                    A[i, j] > 0 && (links[comm[j]] = get(links, comm[j], 0.0) + A[i, j])
                end
                tot[ci] -= k[i]
                best, bestgain = ci, get(links, ci, 0.0) - resolution * tot[ci] * k[i] / twom
                for (c, w) in links
                    gain = w - resolution * tot[c] * k[i] / twom
                    gain > bestgain + 1e-12 && ((best, bestgain) = (c, gain))
                end
                tot[best] += k[i]
                if best != ci
                    comm[i] = best
                    moved = improved = true
                end
            end
        end
        improved || break
        cs = unique(comm)
        idx = Dict(c => t for (t, c) in enumerate(cs))
        newmembers = [Int[] for _ in cs]
        for i in 1:m
            append!(newmembers[idx[comm[i]]], members[i])
        end
        B = zeros(length(cs), length(cs))
        for i in 1:m, j in 1:m
            B[idx[comm[i]], idx[comm[j]]] += A[i, j]
        end
        A, members = B, newmembers
    end
    for (g, mem) in enumerate(members), i in mem
        labels[i] = g
    end
    return labels
end

"""
    modular_groups(W; runs = 10, resolution = 1.0, rng = Random.default_rng())

Groups of nodes of the weighted network `W` that maximise modularity, by the Louvain method,
keeping the best of `runs` randomised runs (the paper ran Gephi's implementation ten times).
Returns `(groups, Q)`, with groups numbered 1, 2, … in order of first appearance.
"""
function modular_groups(W::AbstractMatrix; runs::Integer = 10, resolution::Real = 1.0,
                        rng::AbstractRNG = Random.default_rng())
    best, bestQ = Int[], -Inf
    for _ in 1:runs
        g = louvain_once(rng, W, resolution)
        Q = modularity(W, g; resolution = resolution)
        Q > bestQ + 1e-12 && ((best, bestQ) = (g, Q))
    end
    relabel = Dict{Int,Int}()
    for g in best
        haskey(relabel, g) || (relabel[g] = length(relabel) + 1)
    end
    return [relabel[g] for g in best], bestQ
end

"""
    cophenetic(h::Hclust)

The matrix of cophenetic distances of a hierarchical clustering: for each pair of leaves, the
height at which they are first merged.
"""
function cophenetic(h::Hclust)
    n = length(h.order)
    D = zeros(n, n)
    clusters = Dict{Int,Vector{Int}}()
    members(x) = x < 0 ? [-x] : clusters[x]
    for (s, (a, b)) in enumerate(eachrow(h.merges))
        ma, mb = members(a), members(b)
        for i in ma, j in mb
            D[i, j] = D[j, i] = h.heights[s]
        end
        clusters[s] = vcat(ma, mb)
    end
    return D
end

upper(D) = [D[i, j] for j in axes(D, 2) for i in axes(D, 1) if i < j]

correlation_distance(X) = 1 .- cor(X')   # between rows

"""
    metric_dendrogram(R; distance = :profiles)

Average-linkage hierarchical clustering of the metrics from their mean correlation matrix `R`.
`distance = :profiles` follows the authors' code: rows of `|R|` are treated as correlation
profiles and compared by correlation distance (SciPy `metric='correlation'`).
`distance = :one_minus_rho` follows the paper's text, `d = 1 − ρ`. Returns the `Hclust` and
the cophenetic correlation against both the clustering distance and, as the authors' code
computes it, the Euclidean distance between `|R|` profiles: `(tree, c_clustering, c_code)`.
"""
function metric_dendrogram(R::AbstractMatrix; distance::Symbol = :profiles)
    P = abs.(R)
    D = if distance === :profiles
        correlation_distance(P)
    elseif distance === :one_minus_rho
        1 .- R
    else
        throw(ArgumentError("distance must be :profiles or :one_minus_rho, got $distance"))
    end
    D = (D + D') / 2
    D[diagind(D)] .= 0
    h = hclust(D; linkage = :average)
    C = upper(cophenetic(h))
    E = [norm(P[i, :] - P[j, :]) for j in axes(P, 1) for i in axes(P, 1) if i < j]
    return (tree = h, c_clustering = cor(C, upper(D)), c_code = cor(C, E))
end

"""
    explained_variance(Σ)

For each variable of the covariance matrix `Σ`, the share of the total variance it explains
through its covariances, `Σⱼ Σⱼᵢ² / (tr(Σ) Σᵢᵢ)` (the authors' `get_exp_var_ind`), and its own
share `Σᵢᵢ / tr(Σ)`, as `(exp_var, simple)`.
"""
function explained_variance(Σ::AbstractMatrix)
    T = tr(Σ)
    d = diag(Σ)
    return (vec(sum(abs2, Σ; dims = 1)) ./ (T .* d), d ./ T)
end

unit_ball_volume(d::Integer) = d == 0 ? 1.0 : d == 1 ? 2.0 : 2π / d * unit_ball_volume(d - 2)

"""
    ellipsoid_volume(Σ)

The volume of the covariance ellipsoid of `Σ` after scaling its eigenvalues to sum to one, as
in the authors' `calculate_HV_covellipsoid`: the unit-ball volume times the product of the
square roots of the scaled eigenvalues. Isotropic covariance gives the largest value
(0.8061 in three dimensions).
"""
function ellipsoid_volume(Σ::AbstractMatrix)
    λ = eigvals(Symmetric(Matrix(Σ)))
    λ = max.(λ, 0) ./ sum(λ)
    return unit_ball_volume(size(Σ, 1)) * prod(sqrt, λ)
end

"""
    correlation_types(sc::SizeCorrelations; threshold = 0.1)

The classification of each metric pair by how its correlation changes with community size
(SI §1, Fig. S1), with the authors' size windows: small 5–14 and 5–15, medium 50–59 (the code's
positional slice; 45–55 for the sign test), large 85–95. Returns a `DataFrame` with columns
`a`, `b` and `type`, the panel letter `"A"` to `"G"` (G: unchanged), or `"-"` if unassigned.
"""
function correlation_types(sc::SizeCorrelations; threshold::Real = 0.1)
    m(lo, hi) = dropdims(mean(sc.rho[:, :, sizeindex(sc, lo, hi)]; dims = 3); dims = 3)
    s(lo, hi) = dropdims(std(sc.rho[:, :, sizeindex(sc, lo, hi)]; dims = 3); dims = 3)
    small10, small11, large, med = m(5, 14), m(5, 15), m(85, 95), m(50, 59)
    rows = NamedTuple{(:a, :b, :type),Tuple{String,String,String}}[]
    k = length(sc.labels)
    for i in 1:k, j in (i + 1):k
        cte = abs(small10[i, j] - large[i, j]) < threshold
        ctemedbig = abs(med[i, j] - large[i, j]) < threshold
        decrease = abs(large[i, j]) < abs(small11[i, j])
        increase = abs(large[i, j]) > abs(small11[i, j])
        sign_change = sign(small11[i, j]) != sign(large[i, j]) && abs(large[i, j]) > 0.1 &&
                      abs(small11[i, j]) > 0.1
        t = if cte
            "G"
        elseif sign_change
            ctemedbig ? "D" : "A"
        elseif increase
            ctemedbig ? "E" : "B"
        elseif decrease
            ctemedbig ? "F" : "C"
        else
            "-"
        end
        push!(rows, (a = sc.labels[i], b = sc.labels[j], type = t))
    end
    return DataFrame(rows)
end
