@testset "S-map Jacobians from time series" begin
    rng = Xoshiro(23)
    # a linear map with known continuous-time Jacobian: theta = 0 recovers it
    J = [-0.5 0.2 0.0; -0.3 -0.4 0.1; 0.0 0.25 -0.6]
    Φ = exp(J)
    T = 4000
    X = zeros(T, 3)
    X[1, :] = [1.0, 1.0, 1.0]
    for t in 2:T
        X[t, :] = 1 .+ Φ * (X[t - 1, :] .- 1) .+ 0.05 .* randn(rng, 3)
    end
    js = smap_jacobians(X; theta = 0.0, times = [T ÷ 2])
    @test maximum(abs, matrix(js[1]) - J) < 0.03
    @test resilience(js[1]) ≈ resilience(J) atol = 0.03

    # a noisy Lotka–Volterra community near equilibrium: local fits recover its Jacobian
    a = AlphaMatrix([-1.0 -0.3 0.0; -0.2 -1.0 -0.4; 0.0 -0.1 -1.0], [1.0, 0.9, 0.8])
    N = equilibrium(a)
    Jtrue = matrix(jacobian(a))
    dt, steps, every = 0.01, 200_000, 50
    x = copy(N)
    du = zeros(3)
    rows = Vector{Float64}[]
    for k in 1:steps
        glv_rhs!(du, x, a, 0.0)
        x .= max.(x .+ du .* dt .+ 0.02 .* x .* sqrt(dt) .* randn(rng, 3), 1e-6)
        k % every == 0 && push!(rows, copy(x))
    end
    Y = reduce(vcat, permutedims.(rows))
    near = argmin([norm(Y[t, :] .- N) for t in 1:(size(Y, 1) - 1)])
    js = smap_jacobians(Y; theta = 1.0, dt = dt * every, times = [near])
    @test maximum(abs, matrix(js[1]) - Jtrue) < 0.15
    @test sign.(matrix(js[1])) == sign.(Jtrue) || count(sign.(matrix(js[1])) .!= sign.(Jtrue)) ≤ 2

    @test invariability([1.0, 2.0, 3.0]) ≈ 2.0
end
