@testset "typed matrices and structural stability" begin
    rng = Xoshiro(17)

    # A = −I, r = 1: N* = 1, feasibility domain = positive orthant (Ω^(1/S) = 1/2)
    a = AlphaMatrix(-Matrix(1.0I, 4, 4), ones(4); meta = MatrixMetadata(taxa = ["a", "b", "c", "d"]))
    @test equilibrium(a) ≈ ones(4)
    @test isfeasible(a)
    fd = feasibility_domain(rng, a; samples = 200_000)
    @test isapprox(fd.Ω, 1 / 16; atol = 4fd.se)
    @test isapprox(fd.Ω_scaled, 0.5; atol = 0.01)
    @test feasibility_angle(a) ≈ 0 atol = 1e-6
    @test taxa(a) == ["a", "b", "c", "d"]
    @test_throws ArgumentError AlphaMatrix(ones(2, 3), ones(2))
    @test_throws ArgumentError AlphaMatrix(ones(2, 2), ones(2); meta = MatrixMetadata(taxa = ["x"]))

    # a random competitive community: the Jacobian matches finite differences of the gLV field
    A = -Matrix(1.0I, 5, 5) .- 0.1 .* rand(rng, 5, 5)
    b = AlphaMatrix(A, 1 .+ rand(rng, 5))
    @test isfeasible(b)
    J = jacobian(b)
    N = equilibrium(b)
    f1, f2, Jfd = zeros(5), zeros(5), zeros(5, 5)
    for l in 1:5
        Np, Nm = copy(N), copy(N)
        Np[l] += 1e-6
        Nm[l] -= 1e-6
        glv_rhs!(f1, Np, b, 0.0)
        glv_rhs!(f2, Nm, b, 0.0)
        Jfd[:, l] = (f1 - f2) / 2e-6
    end
    @test matrix(J) ≈ Jfd atol = 1e-7
    @test matrix(per_capita(J)) ≈ A
    @test occursin("Lotka", matrix_metadata(J).derivation)
    @test resilience(J) > 0
    @test keys(jacobian_metrics(J)) == (:Rinf, :R0, :Is, :Amax, :tmax, :S, :absumS, :MRinf, :MR0, :MAmax, :Mtmax)

    # fluxes of a resource–consumer pair reproduce the Lotka–Volterra Jacobian
    r, s, at, e, m = 1.0, 0.5, 2.0, 0.6, 0.3
    R = m / (e * at)
    C = (r - s * R) / at
    F = [0.0 0.0; at*R*C 0.0]
    Jf = jacobian(FluxMatrix(F, [R, C]); efficiency = e, self_regulation = [s * R, 0.0])
    @test matrix(Jf) ≈ [-s*R -at*R; e*at*C 0.0]

    # a self-regulated prey and its predator: sign-stable, press responses fully determined
    sd = SignedDigraph([-1 -1; 1 0])
    @test quasi_sign_stability(rng, sd; samples = 2000).fraction == 1
    pr = press_sign_determinacy(rng, sd; samples = 2000)
    @test all(pr.determinacy .== 1)
    @test pr.sign[1, 1] == 0 && pr.sign[2, 1] == 1     # pressing the prey feeds the predator only
    @test SignedDigraph(J).S == Int8.(sign.(matrix(J)))

    # the typed Jacobian of a simulated community gives the same metrics as the raw functions
    c = nothing
    while c === nothing
        c, _ = assemble_community(rng, 10)
    end
    j = Jacobian(c)
    @test resilience(j) == resilience(jacobian(c.B, c.model))
    @test jacobian_metrics(j).MR0 ≈ median_reactivity(jacobian(c.B, c.model), pulse_covariance(biomass(c)))

    p = ProbabilisticNetwork(fill(0.5, 3, 3))
    @test size(rand(rng, p)) == (3, 3)
end

@testset "Lotka–Volterra assembly and forced trajectories" begin
    # strong competition: species 1 excludes species 2; species 3 is outside the pool
    a = AlphaMatrix([-1.0 -2.0 0.0; -2.0 -1.0 0.0; 0.0 0.0 -1.0], [1.0, 0.8, 1.0])
    ss = steady_state(a, [0.5, 0.5, 0.0])
    @test ss.survivors == [1] && ss.N ≈ [1.0, 0.0, 0.0] && ss.stable && ss.uninvadable
    @test invasion_rates(a, ss.N)[2] ≈ 0.8 - 2.0
    # species 3 can invade but is not in the pool unless named
    @test steady_state(a, [0.5, 0.5, 0.0]; pool = 1:3).survivors == [1, 3]
    @test matrix(subcommunity(a, [1, 3])) == [-1.0 0.0; 0.0 -1.0]

    # a pulse that halves growth for one time unit, then recovery to N* = 1
    b = AlphaMatrix([-1.0;;], [1.0])
    X = glv_trajectory(b, [1.0], [0.0, 1.0, 20.0]; forcing = t -> t < 1 ? 1.0 : 0.0, susceptibility = [-0.5],
                       tstops = [1.0])
    @test X[2, 1] ≈ 0.5 / (1 - 0.5exp(-0.5)) rtol = 1e-6   # logistic with r = 0.5, K = 0.5 from N = 1
    @test X[3, 1] ≈ 1.0 atol = 1e-6
end

@testset "Trophic species of an observed web" begin
    # plant 1; herbivores 2 and 3 both eat only 1; 4 eats both; detritus 5 eaten by 2 and 3
    A = falses(5, 5)
    A[2, 1] = A[3, 1] = A[4, 2] = A[4, 3] = A[2, 5] = A[3, 5] = true
    w = ObservedWeb("toy", "lakes", FoodWeb(A), ["p", "h1", "h2", "c", "d"], [1.0, 4.0, 16.0, 100.0, NaN],
                    BitVector([1, 1, 1, 1, 0]), ["photo", "invertebrate", "invertebrate", "ectotherm", "detritus"])
    t = trophic_species(w)
    @test richness(t.web) == 4
    @test "h1 | h2" in t.taxa
    @test t.mass[findfirst(==("h1 | h2"), t.taxa)] ≈ 8.0
    @test count(t.web.A) == 3                    # H→p, H→d, c→H
    @test trophic_species(t).taxa == t.taxa      # idempotent
end
