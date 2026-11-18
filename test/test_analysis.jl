# Regression against the authors' Python analysis (lab/code/python-reference) on Dataset S2.
@testset "analysis against the Python reference" begin
    df = read_metrics(joinpath(pkgdir(DominguezGarcia2019), "data", "dominguez-garcia2019", "Dataset_S2.csv"))
    @test size(df) == (9380, 29)
    sc = size_correlations(df)
    @test sc.sizes == collect(5:99)
    li = Dict(l => i for (i, l) in enumerate(sc.labels))

    # value covariances do not depend on rank ties, so they match to rounding
    ref = CSV.read(joinpath(REF, "DB_covP.csv"), DataFrame)
    d = maximum(abs(sc.valuecov[li[split(n, ".")[1]], li[split(n, ".")[2]], r] - ref[r, n])
                for r in 1:nrow(ref) for n in names(ref)[2:end])
    @test d < 1e-12

    for (cls, tag) in ("small" => "rSMA", "large" => "rBIG", "medium" => "rMED")
        lo, hi = Dict(SIZE_CLASSES)[cls]
        Σ, _ = class_average(sc, lo, hi; field = :valuecov)
        ev, _ = explained_variance(Σ)
        rk = CSV.read(joinpath(REF, "rank_expvar_full_$tag.csv"), DataFrame)
        @test maximum(abs(ev[li[m]] - v) for (m, v) in zip(rk.metric, rk.exp_var)) < 1e-12
    end

    # the small class has no tie-sensitive differences: the printed cophenetic coefficient
    R, _ = class_average(sc, 5, 15)
    @test round(metric_dendrogram(R).c_code; digits = 6) == 0.871918

    # Figure S1 classification of the 351 pairs (SI §1)
    counts = Dict(k => count(==(k), correlation_types(sc).type) for k in ["A", "B", "C", "D", "E", "F", "G"])
    @test counts == Dict("A" => 7, "B" => 8, "C" => 7, "D" => 5, "E" => 64, "F" => 105, "G" => 155)

    # modularity: Louvain reaches the paper's Q = 0.177 on the medium class
    Rm, _ = class_average(sc, 45, 55)
    g, Q = modular_groups(metric_network(Rm); runs = 50, rng = Xoshiro(1))
    @test round(Q; digits = 3) == 0.177
    @test modularity(metric_network(Rm), g) ≈ Q

    @test ellipsoid_volume(Matrix(1.0I, 3, 3)) ≈ 0.8061 atol = 1e-4
end
