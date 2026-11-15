@testset "Jacobian-based metrics" begin
    # J = −a I: not reactive, resilience a, Ĵ = −2a I so Is = a, sensitivity matrix I/a
    a = 0.7
    J = -a * Matrix(1.0I, 4, 4)
    @test reactivity(J) ≈ -a
    @test resilience(J) ≈ a
    @test max_amplification(J) == (1.0, 0.0)
    @test stochastic_invariability(J) ≈ a
    @test press_sensitivity(J) ≈ 4 / a
    @test total_sensitivity(J) ≈ 4 / a

    # a non-normal, stable but reactive matrix
    Jn = [-1.0 10.0; 0.0 -2.0]
    @test reactivity(Jn) > 0
    Amax, tmax = max_amplification(Jn)
    @test Amax > 1 && tmax > 0
    @test Amax ≥ amplification_envelope(Jn, tmax) - 1e-12
    @test isapprox(round(tmax; digits = 1), tmax; atol = 1e-9)
    # the ∞-norm envelope (published Amax): maximum absolute row sum of exp(Jt), on the same grid
    Ainf, tinf = max_amplification(Jn; p = Inf)
    grid = 0.1:0.1:20
    @test Ainf ≈ maximum(t -> maximum(sum(abs, exp(Jn * t); dims = 2)), grid)
    @test amplification_envelope(Jn, tinf; p = Inf) ≈ Ainf
    @test max_amplification(J; p = Inf) == (1.0, 0.0)

    # dense and iterative invariability agree
    rng = Xoshiro(3)
    for n in (8, 25)
        R = randn(rng, n, n) - 2n * I
        @test isapprox(stochastic_invariability(R; dense_max = 100), stochastic_invariability(R; dense_max = 0);
                       rtol = 1e-6)
    end

    # median metrics with biomass-weighted pulses
    B = [1.0, 0.5]
    C = pulse_covariance(B)
    @test median_reactivity(Jn, C) ≈ (-1.0 * 1 + -2.0 * 0.25) / 1.25
    @test median_displacement(Jn, C, 0.0) ≈ 1
    @test median_resilience(-a * Matrix(1.0I, 2, 2), C) ≈ a rtol = 1e-9
end

@testset "Median amplification after an initial decline" begin
    J = [-1.0 0.0; 5.0 -0.1]            # non-normal: species 1 feeds a slow species 2
    C = Matrix(1.0I, 2, 2)
    @test median_reactivity(J, C) < 0   # the median displacement falls at first
    MAmax, Mtmax = median_max_amplification(J, C)
    @test MAmax > 2
    @test 1 < Mtmax < 10
    @test median_max_amplification(-Matrix(1.0I, 2, 2), C) == (1.0, 0.0)
end
