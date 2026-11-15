# Structural stability of Lotka–Volterra coefficient matrices (Rohr, Saavedra & Bascompte 2014;
# Saavedra et al. 2017) and qualitative analysis of sign patterns (Allesina & Pascual 2008;
# Dambacher et al. 2002, 2003). These are the analyses that need only a matrix, and so apply to
# empirical α matrices and sign patterns as well as to simulated communities.

"""
    feasibility_domain([rng], a::AlphaMatrix; samples = 100_000)

The size of the feasibility domain of `a` (Saavedra et al. 2017): the fraction `Ω` of
directions of the growth-rate vector `r` (uniform on the unit sphere) for which the interior
equilibrium `−A⁻¹ r` is positive, estimated from `samples` random directions. Returns
`(Ω, Ω^(1/S), se)`, with `se` the Monte Carlo standard error of `Ω`. `Ω^(1/S)` is comparable
across community sizes; for `A = −I` it is 1/2.
"""
function feasibility_domain(rng::AbstractRNG, a::AlphaMatrix; samples::Integer = 100_000)
    S = nspecies(a)
    F = lu(a.A)
    hits = 0
    batch = 10_000
    done = 0
    while done < samples
        m = min(batch, samples - done)
        X = -(F \ randn(rng, S, m))
        hits += count(k -> all(>(0), view(X, :, k)), 1:m)
        done += m
    end
    Ω = hits / samples
    return (Ω = Ω, Ω_scaled = Ω^(1 / S), se = sqrt(Ω * (1 - Ω) / samples))
end
feasibility_domain(a::AlphaMatrix; kwargs...) = feasibility_domain(Random.default_rng(), a; kwargs...)

"""
    structural_vector(a::AlphaMatrix)

The centroid of the feasibility cone of `a`: the normalised sum of the unit columns of `−A`
(Saavedra et al. 2017). Growth-rate vectors along it are the most robustly feasible.
"""
function structural_vector(a::AlphaMatrix)
    V = -a.A
    v = vec(sum(V ./ sqrt.(sum(abs2, V; dims = 1)); dims = 2))
    return v ./ norm(v)
end

"""
    feasibility_angle(a::AlphaMatrix)

The angle in degrees between the growth-rate vector `a.r` and the [`structural_vector`](@ref)
of `a`: a measure of how far the community is from the centre of its feasibility domain
(larger means closer to losing feasibility).
"""
function feasibility_angle(a::AlphaMatrix)
    c = clamp(dot(a.r, structural_vector(a)) / norm(a.r), -1, 1)
    return rad2deg(acos(c))
end

function random_from_signs(rng::AbstractRNG, S::AbstractMatrix)
    return S .* rand(rng, size(S)...)
end

"""
    quasi_sign_stability([rng], s::SignedDigraph; samples = 10_000)

The fraction of matrices with sign pattern `s` and magnitudes uniform on (0, 1) that are stable
(every eigenvalue with negative real part): quasi-sign-stability (Allesina & Pascual 2008).
A diagonal of zeros means no self-regulation. Returns `(fraction, se)`.
"""
function quasi_sign_stability(rng::AbstractRNG, s::SignedDigraph; samples::Integer = 10_000)
    stable = count(_ -> maximum(real, eigvals(random_from_signs(rng, s.S))) < 0, 1:samples)
    p = stable / samples
    return (fraction = p, se = sqrt(p * (1 - p) / samples))
end
quasi_sign_stability(s::SignedDigraph; kwargs...) = quasi_sign_stability(Random.default_rng(), s; kwargs...)

"""
    press_sign_determinacy([rng], s::SignedDigraph; samples = 10_000)

Qualitative press-perturbation analysis (Dambacher et al. 2003; the simulation approach of
QPress): draw matrices with sign pattern `s` and magnitudes uniform on (0, 1), keep the stable
ones, and for each compute the equilibrium response to a sustained increase in each species'
growth, `−J⁻¹`. Returns `(sign, determinacy, stable)`: the most common sign of each response
(row: responding species, column: pressed species; 0 when the response is structurally zero,
below 1e-10 of the largest), the fraction of stable draws with that outcome, and the fraction
of draws that were stable.
"""
function press_sign_determinacy(rng::AbstractRNG, s::SignedDigraph; samples::Integer = 10_000)
    n = nspecies(s)
    pos, neg, kept = zeros(Int, n, n), zeros(Int, n, n), 0
    for _ in 1:samples
        J = random_from_signs(rng, s.S)
        maximum(real, eigvals(J)) < 0 || continue
        R = -inv(J)
        tol = 1e-10 * maximum(abs, R)
        pos .+= R .> tol
        neg .+= R .< -tol
        kept += 1
    end
    kept == 0 && return (sign = zeros(Int8, n, n), determinacy = fill(NaN, n, n), stable = 0.0)
    zero_ = kept .- pos .- neg
    sgn = [zero_[i] ≥ max(pos[i], neg[i]) ? Int8(0) : Int8(sign(pos[i] - neg[i])) for i in eachindex(pos)]
    det = max.(pos, neg, zero_) ./ kept
    return (sign = sgn, determinacy = det, stable = kept / samples)
end
press_sign_determinacy(s::SignedDigraph; kwargs...) = press_sign_determinacy(Random.default_rng(), s; kwargs...)
