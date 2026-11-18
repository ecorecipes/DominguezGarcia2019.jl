@testset "published Amax convention" begin
    # a stable, reactive, non-normal Jacobian: the ∞-norm envelope on the grid
    Jn = [-1.0 10.0; 0.0 -2.0]
    op = ContinuousJacobian(Jn)
    Ainf, tinf = published_max_amplification(op)
    @test Ainf ≈ maximum(t -> maximum(sum(abs, exp(Jn * t); dims = 2)), 0.1:0.1:20)
    @test amplification_envelope(op, tinf; p = Inf) ≈ Ainf
    @test isapprox(round(tinf; digits = 1), tinf; atol = 1e-9)
    # not reactive in the 2-norm: Amax = 1 exactly, as in Dataset S1 (Amax = 1 iff R0 ≤ 0)
    @test published_max_amplification(ContinuousJacobian(-0.7 * Matrix(1.0I, 3, 3))) == (1.0, 0.0)
    # R0 ≤ 0 but the ∞-norm logarithmic norm is positive: the published gate still gives 1,
    # where ResilienceCore's ∞-norm scan finds amplification
    Jg = [-1.0 0.0; 1.8 -1.0]
    @test reactivity(ContinuousJacobian(Jg)) ≤ 0
    @test published_max_amplification(ContinuousJacobian(Jg)) == (1.0, 0.0)
    @test max_amplification(ContinuousJacobian(Jg); p = Inf)[1] > 1
end

@testset "communities of the 10,000-row ensemble are reproduced from their seeds" begin
    # ensemble_v4 was computed with EcologicalResilience's own plain-matrix metrics (Apple
    # silicon, one BLAS thread); ResilienceCore's metrics give the stored Jacobian-based
    # columns, up to the rounding of LAPACK on other platforms and thread counts
    ens = CSV.read(joinpath(pkgdir(DominguezGarcia2019), "data", "simulated", "ensemble_v4.csv"), DataFrame)
    rows = first(sort(ens[ens.iniN .<= 20, :], :seconds), 8)
    for r in eachrow(rows)
        rng = Xoshiro(parse(Int, r.network[4:end]))
        S = rand(rng, 5:110)
        c = nothing
        while c === nothing
            c, _ = assemble_community(rng, S)
        end
        op = community_jacobian(c)
        pub, si = amplification_metrics(c), amplification_metrics(c; convention = :si)
        same(a, b) = isapprox(a, b; rtol = 1e-6)
        @test same(resilience(op), r.Rinf) && same(reactivity(op), r.R0) && same(stochastic_invariability(op), r.Is)
        @test same(press_sensitivity(op), r.L1DBRMt_s) && same(total_sensitivity(op), r.absumS)
        @test all(same(getfield(pub, k), r[k]) for k in keys(pub))
        @test all(same(getfield(si, k), r[Symbol(k, "_si")]) for k in keys(si))
    end
end
