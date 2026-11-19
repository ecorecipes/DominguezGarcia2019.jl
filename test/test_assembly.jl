@testset "community assembly and press experiments" begin
    rng = Xoshiro(2024)
    c, status = nothing, :none
    while c === nothing
        c, status = assemble_community(rng, 12)
    end
    @test status == :ok
    op = community_jacobian(c)
    @test op.J == jacobian(c.B, c.model) && op.equilibrium == biomass(c)
    @test maximum(real, eigvals(op.J)) < 0
    @test all(biomass(c) .> c.model.threshold[c.model.alive])
    @test isempty(isolated_plants(c.model))
    @test richness(c) ≥ floor(Int, 0.9 * c.initial_richness)   # Dataset S1: N ≥ floor(0.9 iniN)

    opts = PerturbationOptions(sequences = 3)
    TM = global_mortality_tolerance(c; options = opts)
    @test isapprox(round((TM - 1.001) / 0.1), (TM - 1.001) / 0.1; atol = 1e-9)
    loc = local_mortality_tolerances(c; options = opts)
    @test length(loc) == richness(c) && all(>(0), loc)
    rem = removal_responses(c; options = opts)
    @test all(r -> r.survivors ≤ richness(c) - 1, rem)

    r = community_record(rng, c; options = opts)
    @test keys(r) == Tuple(RAW_COLUMNS)
    pub, si = amplification_metrics(c), amplification_metrics(c; convention = :si)
    @test all(getfield(r, k) == getfield(pub, k) for k in keys(pub))
    @test si.Amax ≈ max_amplification(op)[1] && pub.Amax ≈ published_max_amplification(op)[1]
    @test si.MR0 ≈ median_reactivity(op)
    @test_throws ArgumentError amplification_metrics(c; convention = :other)
    # published weights: the first N entries of the full biomass vector; aligned when nothing died
    N = richness(c)
    @test published_pulse_covariance(c) == Diagonal(c.B[1:N] .^ 2)
    alive = c.model.alive
    full = Community(FoodWebModel(FoodWeb(c.model.web.A[alive, alive]), c.model.params, c.model.exponents),
                     c.B[alive], N, N)
    @test published_pulse_covariance(full) == pulse_covariance(biomass(full))
    # the published global press: an absolute mortality rise of 0.055 on every species
    @test published_press_factors(c) ≈ 1 .+ 0.055 ./ c.model.base_mortality
    @test r.RM_p == press_response(c, published_press_factors(c); options = opts)[1]
    prepared = prepare_metrics(DataFrame([r]))
    @test names(prepared) == [METRICS; "N"]
    # ResilienceCore's named profile gives the paper's 27 metrics under the SI conventions
    prof = DominguezGarcia2019.ResilienceCore.stability_profile(Xoshiro(7), c; options = opts)
    rsi = prepare_metrics(DataFrame([community_record(Xoshiro(7), c; options = opts, convention = :si)]))
    k = N / (N - 1)
    sgn = [-1, 1, -1, -1, 1, -1, -1, -1, 1, -1, -1, 1, 1, 1, -k, 1, -k, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
    @test length(prof) == length(METRICS) == length(sgn)
    @test all(isapprox(rsi[1, m], s * v; rtol = 1e-10, atol = 1e-14, nans = true)
              for (m, s, v) in zip(METRICS, sgn, values(prof)))
end
