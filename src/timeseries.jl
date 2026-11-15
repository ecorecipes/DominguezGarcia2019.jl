# Interaction matrices from time series (route 2 of the interaction-matrix review): local
# Jacobians by the multivariate S-map (Sugihara 1994; Deyle et al. 2016), optionally
# regularised (Cenci, Sugihara & Saavedra 2019), and stability measures computed from data
# alone. The continuous-time Jacobians feed every Jacobian metric of the package.

"""
    smap_jacobians(X; theta = 2.0, ridge = 0.0, dt = 1.0, times = 1:size(X, 1) - 1,
                   continuous = true, meta = MatrixMetadata())

Local Jacobians from a multivariate time series `X` (rows are times at spacing `dt`, columns
species) by the S-map. For each target time `t*`, every `xᵢ(t + 1)` is regressed linearly on
`x(t)` (with intercept), weighting each past state by `exp(−theta d / d̄)`, where `d` is its
distance from `x(t*)` and `d̄` the mean distance; `ridge > 0` adds an L2 penalty on the slopes.
The slopes form the one-step map Jacobian `Φ(t*)`. With `continuous = true` it is converted to
a continuous-time [`Jacobian`](@ref) by `log(Φ)/dt` (real part); otherwise `Φ` is returned in
a `Jacobian` (map) and only its spectral radius is meaningful for stability. `theta = 0` gives
one global linear (VAR(1)) fit. Returns a vector of `Jacobian`s, one per element of `times`,
each with the state `x(t*)` as its equilibrium field (a reference state, not an equilibrium).
"""
function smap_jacobians(X::AbstractMatrix; theta::Real = 2.0, ridge::Real = 0.0, dt::Real = 1.0,
                        times = 1:(size(X, 1) - 1), continuous::Bool = true,
                        meta::MatrixMetadata = MatrixMetadata())
    T, S = size(X)
    T > S + 2 || throw(ArgumentError("need more time points ($T) than species + 2 ($(S + 2))"))
    lib = 1:(T - 1)
    Xl = Matrix{Float64}(X[lib, :])
    Yl = Matrix{Float64}(X[lib .+ 1, :])
    D = hcat(ones(length(lib)), Xl)
    P = Diagonal([0.0; fill(Float64(ridge), S)])
    note = string("S-map, theta = ", theta, ", ridge = ", ridge, ", ", T, " time points",
                  continuous ? ", continuous-time by matrix logarithm" : ", one-step map")
    out = Jacobian[]
    for t in times
        x = X[t, :]
        d = [norm(view(Xl, k, :) .- x) for k in lib]
        dbar = mean(d)
        w = dbar > 0 ? exp.(-theta .* d ./ dbar) : ones(length(lib))
        Wd = D .* w
        β = (D' * Wd + P) \ (Wd' * Yl)
        Φ = Matrix(β[2:end, :]')
        J = continuous ? real.(log(complex.(Φ))) ./ dt : Φ
        push!(out, Jacobian(J, x; meta = derived(meta, note)))
    end
    return out
end

"""
    invariability(x)

Temporal invariability `mean(x) / std(x)` of a series (Haegeman et al. 2016), e.g. of total
community biomass; the data counterpart of stochastic invariability.
"""
invariability(x::AbstractVector) = mean(x) / std(x)
