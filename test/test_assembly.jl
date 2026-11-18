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
    prepared = prepare_metrics(DataFrame([r]))
    @test names(prepared) == [METRICS; "N"]
end
