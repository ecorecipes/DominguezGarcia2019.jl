# Typed interaction matrices. "Interaction matrix" names several different objects (Terry et
# al. 2025): Lotka–Volterra coefficients, the Jacobian at an equilibrium, the per-capita
# community matrix, flow matrices and sign patterns. Software written for one silently accepts
# another. Here each is its own type, conversions are explicit functions that ask for what they
# need, and every matrix carries a metadata record (lab/wiki/project/interaction-matrices.md).

"""
    MatrixMetadata(; taxa = String[], state_variables, units, sign_convention, spatial_bounds,
                     temporal_bounds, derivation, uncertainty, software, licence)

The minimal metadata record proposed for deposited interaction matrices (interaction-matrix
review, lab/papers/interaction): the species (`taxa`, in row order), what the state variables
are, units, the sign convention, spatial and temporal bounds, how the entries were derived,
how their uncertainty is represented, the software used and the licence. All fields are
strings except `taxa`; unknown fields default to `"unspecified"`.
"""
Base.@kwdef struct MatrixMetadata
    taxa::Vector{String} = String[]
    state_variables::String = "unspecified"
    units::String = "unspecified"
    sign_convention::String = "entry (i, j) is the effect of species j on species i"
    spatial_bounds::String = "unspecified"
    temporal_bounds::String = "unspecified"
    derivation::String = "unspecified"
    uncertainty::String = "point estimate"
    software::String = "EcologicalResilience.jl"
    licence::String = "unspecified"
end

# metadata with `derivation` extended by a note on the conversion just applied
derived(m::MatrixMetadata, note::AbstractString) =
    MatrixMetadata(m.taxa, m.state_variables, m.units, m.sign_convention, m.spatial_bounds, m.temporal_bounds,
                   m.derivation == "unspecified" ? note : string(m.derivation, "; ", note), m.uncertainty,
                   m.software, m.licence)

"""
    InteractionMatrix

Supertype of the typed interaction matrices: [`AlphaMatrix`](@ref), [`Jacobian`](@ref),
[`CommunityMatrix`](@ref), [`FluxMatrix`](@ref), [`SignedDigraph`](@ref) and
[`ProbabilisticNetwork`](@ref). Use [`matrix`](@ref) for the numbers, [`matrix_metadata`](@ref) and
[`taxa`](@ref) for the record, and `size`/`nspecies` for the dimension.
"""
abstract type InteractionMatrix end

"""
    matrix(x::InteractionMatrix)

The numeric matrix held by `x` (rows are affected species, columns are acting species).
"""
function matrix end

"""
    matrix_metadata(x::InteractionMatrix)

The [`MatrixMetadata`](@ref) record of `x`.
"""
matrix_metadata(x::InteractionMatrix) = x.meta

"""
    taxa(x::InteractionMatrix)

Species names in row order (`"sp1"`, `"sp2"`, … if none were given).
"""
taxa(x::InteractionMatrix) = isempty(x.meta.taxa) ? ["sp$i" for i in 1:nspecies(x)] : x.meta.taxa

"""
    nspecies(x::InteractionMatrix)

The number of species (rows) of `x`.
"""
nspecies(x::InteractionMatrix) = size(matrix(x), 1)
Base.size(x::InteractionMatrix) = size(matrix(x))

function checksquare(M, name)
    size(M, 1) == size(M, 2) || throw(ArgumentError("$name must be square, got $(size(M))"))
    return size(M, 1)
end
function checkmeta(meta::MatrixMetadata, n)
    isempty(meta.taxa) || length(meta.taxa) == n ||
        throw(ArgumentError("metadata lists $(length(meta.taxa)) taxa for $n species"))
    return meta
end

"""
    AlphaMatrix(A, r; meta = MatrixMetadata())

Generalised Lotka–Volterra coefficients: `dNᵢ/dt = Nᵢ (rᵢ + Σⱼ Aᵢⱼ Nⱼ)`, with per-capita effects
`A` (negative diagonal for self-regulation) and intrinsic growth rates `r`. Other
parameterisations (e.g. annual-plant `λᵢ/(1 + Σⱼ αᵢⱼ Nⱼ)`) must be converted first; the
conversion belongs in `meta.derivation`.
"""
struct AlphaMatrix <: InteractionMatrix
    A::Matrix{Float64}
    r::Vector{Float64}
    meta::MatrixMetadata
    function AlphaMatrix(A::AbstractMatrix, r::AbstractVector; meta::MatrixMetadata = MatrixMetadata())
        n = checksquare(A, "A")
        length(r) == n || throw(ArgumentError("r has $(length(r)) entries for $n species"))
        return new(Matrix{Float64}(A), Vector{Float64}(r), checkmeta(meta, n))
    end
end
matrix(a::AlphaMatrix) = a.A

"""
    glv_rhs!(du, N, a::AlphaMatrix, t)

The generalised Lotka–Volterra vector field of `a` in the in-place SciML form.
"""
function glv_rhs!(du, N, a::AlphaMatrix, t)
    mul!(du, a.A, N)
    @. du = N * (a.r + du)
    return du
end

"""
    equilibrium(a::AlphaMatrix)

The interior equilibrium `N* = −A⁻¹ r` (all species present), which may have negative entries.
"""
equilibrium(a::AlphaMatrix) = -(a.A \ a.r)

"""
    isfeasible(a::AlphaMatrix)

Whether the interior equilibrium of `a` has every species at positive abundance.
"""
isfeasible(a::AlphaMatrix) = all(>(0), equilibrium(a))

"""
    Jacobian(J, equilibrium = Float64[]; meta = MatrixMetadata())
    Jacobian(c::Community)

The Jacobian (population-level community matrix) of a dynamical system at an equilibrium, with
that equilibrium when known. All Jacobian-based metrics accept it; the median metrics then use
pulses proportional to `equilibrium`.
"""
struct Jacobian <: InteractionMatrix
    J::Matrix{Float64}
    equilibrium::Vector{Float64}
    meta::MatrixMetadata
    function Jacobian(J::AbstractMatrix, equilibrium::AbstractVector = Float64[]; meta::MatrixMetadata = MatrixMetadata())
        n = checksquare(J, "J")
        isempty(equilibrium) || length(equilibrium) == n ||
            throw(ArgumentError("equilibrium has $(length(equilibrium)) entries for $n species"))
        return new(Matrix{Float64}(J), Vector{Float64}(equilibrium), checkmeta(meta, n))
    end
end
matrix(j::Jacobian) = j.J

function Jacobian(c::Community)
    idx = findall(c.model.alive)
    meta = MatrixMetadata(; taxa = ["sp$i" for i in idx], state_variables = "biomass density",
                          units = "per unit time", derivation = "analytic Jacobian of the bioenergetic model at a polished equilibrium")
    return Jacobian(jacobian(c.B, c.model), biomass(c); meta = meta)
end

"""
    jacobian(a::AlphaMatrix; N = equilibrium(a))

The Jacobian of the Lotka–Volterra system `a` at abundances `N`: `diag(N) A` at an interior
equilibrium (for other `N` the general form `diag(r + A N) + diag(N) A` is used).
"""
function jacobian(a::AlphaMatrix; N::AbstractVector = equilibrium(a))
    J = Diagonal(N) * a.A + Diagonal(a.r .+ a.A * N)
    return Jacobian(J, N; meta = derived(a.meta, "Jacobian of generalised Lotka–Volterra at N"))
end

"""
    CommunityMatrix(M, equilibrium; meta = MatrixMetadata())

Per-capita effects at an equilibrium: `Mᵢⱼ = ∂(Ṅᵢ/Nᵢ)/∂Nⱼ`. For Lotka–Volterra it is the
coefficient matrix `A`. Convert with `jacobian(m) = diag(N*) M` and `per_capita(J)`.
"""
struct CommunityMatrix <: InteractionMatrix
    M::Matrix{Float64}
    equilibrium::Vector{Float64}
    meta::MatrixMetadata
    function CommunityMatrix(M::AbstractMatrix, equilibrium::AbstractVector; meta::MatrixMetadata = MatrixMetadata())
        n = checksquare(M, "M")
        length(equilibrium) == n || throw(ArgumentError("equilibrium has $(length(equilibrium)) entries for $n species"))
        return new(Matrix{Float64}(M), Vector{Float64}(equilibrium), checkmeta(meta, n))
    end
end
matrix(m::CommunityMatrix) = m.M

jacobian(m::CommunityMatrix) =
    Jacobian(Diagonal(m.equilibrium) * m.M, m.equilibrium; meta = derived(m.meta, "Jacobian = diag(N*) M"))

"""
    per_capita(j::Jacobian)

The per-capita [`CommunityMatrix`](@ref) `diag(N*)⁻¹ J`; needs the equilibrium, and every
species present.
"""
function per_capita(j::Jacobian)
    isempty(j.equilibrium) && throw(ArgumentError("per_capita needs the equilibrium of the Jacobian"))
    all(>(0), j.equilibrium) || throw(ArgumentError("per_capita needs every species at positive abundance"))
    return CommunityMatrix(Diagonal(1 ./ j.equilibrium) * j.J, j.equilibrium; meta = derived(j.meta, "M = diag(N*)⁻¹ J"))
end

"""
    FluxMatrix(F, biomass; meta = MatrixMetadata())

Feeding fluxes at a steady state: `F[i, j]` is the rate at which consumer `i` removes biomass
of resource `j` (e.g. from Ecopath or fluxweb), with the standing `biomass` of every species.
[`jacobian`](@ref) turns it into a Jacobian under stated assumptions.
"""
struct FluxMatrix <: InteractionMatrix
    F::Matrix{Float64}
    biomass::Vector{Float64}
    meta::MatrixMetadata
    function FluxMatrix(F::AbstractMatrix, biomass::AbstractVector; meta::MatrixMetadata = MatrixMetadata())
        n = checksquare(F, "F")
        length(biomass) == n || throw(ArgumentError("biomass has $(length(biomass)) entries for $n species"))
        all(≥(0), F) || throw(ArgumentError("fluxes must be nonnegative"))
        return new(Matrix{Float64}(F), Vector{Float64}(biomass), checkmeta(meta, n))
    end
end
matrix(f::FluxMatrix) = f.F

"""
    jacobian(f::FluxMatrix; efficiency, self_regulation)

The Jacobian implied by steady-state fluxes when feeding is mass action, `Fᵢⱼ = aᵢⱼ Bᵢ Bⱼ`
(de Ruiter et al. 1995): the effect of resource `j` on consumer `i` is `eⱼ Fᵢⱼ / Bⱼ`, the effect
of consumer `j` on resource `i` is `−Fⱼᵢ / Bⱼ`, and the diagonal is `−self_regulation` (which
fluxes cannot determine). `efficiency` is a vector over resources (or a scalar);
`self_regulation` a vector (or scalar) of nonnegative rates.
"""
function jacobian(f::FluxMatrix; efficiency, self_regulation)
    n = nspecies(f)
    e = efficiency isa Number ? fill(Float64(efficiency), n) : collect(Float64, efficiency)
    s = self_regulation isa Number ? fill(Float64(self_regulation), n) : collect(Float64, self_regulation)
    J = zeros(n, n)
    for i in 1:n, j in 1:n
        i == j && continue
        J[i, j] = e[j] * f.F[i, j] / f.biomass[j] - f.F[j, i] / f.biomass[j]
    end
    J[diagind(J)] .= -s
    note = "Jacobian from mass-action fluxes (de Ruiter et al. 1995), efficiency and self-regulation supplied"
    return Jacobian(J, f.biomass; meta = derived(f.meta, note))
end

"""
    SignedDigraph(S; meta = MatrixMetadata())
    SignedDigraph(x::InteractionMatrix)

A sign pattern (entries −1, 0, +1) of the effects of species on each other, as used in loop
analysis and qualitative press-perturbation analysis. The second form takes the signs of the
Jacobian or community matrix of `x`.
"""
struct SignedDigraph <: InteractionMatrix
    S::Matrix{Int8}
    meta::MatrixMetadata
    function SignedDigraph(S::AbstractMatrix; meta::MatrixMetadata = MatrixMetadata())
        n = checksquare(S, "S")
        all(in((-1, 0, 1)), S) || throw(ArgumentError("a sign pattern has entries −1, 0 or +1"))
        return new(Matrix{Int8}(S), checkmeta(meta, n))
    end
end
matrix(s::SignedDigraph) = s.S
SignedDigraph(x::Union{Jacobian,CommunityMatrix,AlphaMatrix}) =
    SignedDigraph(Int8.(sign.(matrix(x))); meta = derived(x.meta, "signs of $(nameof(typeof(x)))"))

"""
    ProbabilisticNetwork(P; meta = MatrixMetadata())

Probabilities `P[i, j]` that species `j` affects species `i` (for example, that `i` eats `j`),
e.g. from sampling-corrected or predicted links. `rand(rng, p)` draws a `BitMatrix` of links.
"""
struct ProbabilisticNetwork <: InteractionMatrix
    P::Matrix{Float64}
    meta::MatrixMetadata
    function ProbabilisticNetwork(P::AbstractMatrix; meta::MatrixMetadata = MatrixMetadata())
        n = checksquare(P, "P")
        all(p -> 0 ≤ p ≤ 1, P) || throw(ArgumentError("link probabilities must lie in [0, 1]"))
        return new(Matrix{Float64}(P), checkmeta(meta, n))
    end
end
matrix(p::ProbabilisticNetwork) = p.P
Base.rand(rng::AbstractRNG, p::ProbabilisticNetwork) = rand(rng, size(p.P)...) .< p.P

# Jacobian metrics on the typed Jacobian
for f in (:reactivity, :resilience, :max_amplification, :stochastic_invariability, :sensitivity_matrix,
          :press_sensitivity, :total_sensitivity)
    @eval $f(j::Jacobian; kwargs...) = $f(j.J; kwargs...)
end
amplification_envelope(j::Jacobian, t::Real) = amplification_envelope(j.J, t)
function pulses(j::Jacobian)
    isempty(j.equilibrium) && throw(ArgumentError("median metrics need the equilibrium of the Jacobian"))
    return pulse_covariance(j.equilibrium)
end
median_reactivity(j::Jacobian) = median_reactivity(j.J, pulses(j))
median_max_amplification(j::Jacobian; kwargs...) = median_max_amplification(j.J, pulses(j); kwargs...)
median_resilience(j::Jacobian; kwargs...) = median_resilience(j.J, pulses(j); kwargs...)

"""
    jacobian_metrics(j::Jacobian)

The Jacobian-based stability metrics of Domínguez-García et al. (2019) for any Jacobian, in the
raw sign convention of their dataset: `Rinf`, `R0`, `Is`, `Amax`, `tmax`, `S` (community
sensitivity) and `absumS`, and, when the equilibrium is known, `MRinf`, `MR0`, `MAmax`,
`Mtmax`.
"""
function jacobian_metrics(j::Jacobian)
    Amax, tmax = max_amplification(j)
    base = (Rinf = resilience(j), R0 = reactivity(j), Is = stochastic_invariability(j), Amax = Amax, tmax = tmax,
            S = press_sensitivity(j), absumS = total_sensitivity(j))
    isempty(j.equilibrium) && return base
    MAmax, Mtmax = median_max_amplification(j)
    return merge(base, (MRinf = median_resilience(j), MR0 = median_reactivity(j), MAmax = MAmax, Mtmax = Mtmax))
end

"""
    subcommunity(a::AlphaMatrix, idx)

The Lotka–Volterra system of the species `idx` of `a` alone.
"""
function subcommunity(a::AlphaMatrix, idx)
    t = isempty(a.meta.taxa) ? String[] : a.meta.taxa[idx]
    m = a.meta
    meta = MatrixMetadata(t, m.state_variables, m.units, m.sign_convention, m.spatial_bounds, m.temporal_bounds,
                          m.derivation, m.uncertainty, m.software, m.licence)
    return AlphaMatrix(a.A[idx, idx], a.r[idx]; meta = meta)
end

"""
    invasion_rates(a::AlphaMatrix, N)

Per-capita growth rates `r + A N` of every species at abundances `N`; for species absent from
`N` these are their invasion growth rates.
"""
invasion_rates(a::AlphaMatrix, N::AbstractVector) = a.r .+ a.A * N

"""
    steady_state(a::AlphaMatrix, N0; pool = findall(>(0), N0), tmax = 1e4, threshold = 1e-8,
                 reseed = 1e-6, rounds = 20)

The assembled state of the Lotka–Volterra system `a` from `N0`: integrate to `tmax`, remove
species below `threshold`, and solve for the survivors' interior equilibrium; while some species
of the regional `pool` (by default those present in `N0`) could invade (positive invasion
rate), reintroduce them at density `reseed` and repeat (up to `rounds` times). Returns
`(N, survivors, stable, uninvadable)`: the equilibrium (zeros for absent species), the
survivors' indices, whether their Jacobian is stable, and whether every absent species of the
pool has a negative invasion rate.
"""
function steady_state(a::AlphaMatrix, N0::AbstractVector; pool = findall(>(0), N0), tmax::Real = 1e4,
                      threshold::Real = 1e-8, reseed::Real = 1e-6, rounds::Integer = 20)
    x0 = collect(Float64, N0)
    local N, surv, stable, uninvadable
    for _ in 1:rounds
        sol = SciMLBase.solve(ODEProblem(glv_rhs!, x0, (0.0, Float64(tmax)), a), Tsit5(); abstol = 1e-12,
                              reltol = 1e-10, save_everystep = false)
        x = sol.u[end]
        surv = findall(>(threshold), x)
        N = zeros(nspecies(a))
        if !isempty(surv)
            sub = subcommunity(a, surv)
            N[surv] = equilibrium(sub)
        end
        stable = all(>(0), N[surv]) && (isempty(surv) || resilience(jacobian(subcommunity(a, surv))) > 0)
        absent = setdiff(pool, surv)
        g = invasion_rates(a, N)[absent]
        uninvadable = all(<(0), g)
        (uninvadable || !stable) && break
        x0 = copy(N)
        x0[absent[g .> 0]] .= reseed
    end
    return (N = N, survivors = surv, stable = stable, uninvadable = uninvadable)
end

"""
    glv_trajectory(a::AlphaMatrix, N0, times; forcing = nothing, susceptibility = nothing,
                   tstops = Float64[])

Abundances of the Lotka–Volterra system `a` at `times` from `N0` at `times[1]`, optionally
with a perturbation that adds `susceptibility[i] × forcing(t)` to species i's per-capita
growth rate (the antibiotic model of Stein et al. 2013); give the forcing's switching times
as `tstops`. Returns a matrix with one row per time.
"""
function glv_trajectory(a::AlphaMatrix, N0::AbstractVector, times::AbstractVector; forcing = nothing,
                        susceptibility = nothing, tstops = Float64[])
    f! = if forcing === nothing
        glv_rhs!
    else
        ε = collect(Float64, susceptibility)
        (du, N, p, t) -> (glv_rhs!(du, N, p, t); du .+= N .* ε .* forcing(t); du)
    end
    prob = ODEProblem(f!, collect(Float64, N0), (Float64(first(times)), Float64(last(times))), a)
    sol = SciMLBase.solve(prob, Tsit5(); abstol = 1e-12, reltol = 1e-10, saveat = times, tstops = tstops)
    return permutedims(reduce(hcat, sol.u))
end
