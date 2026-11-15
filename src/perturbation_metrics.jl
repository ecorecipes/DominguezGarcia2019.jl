# Stability metrics from simulated press perturbations (SI §10): tolerance to increased
# mortality, resistance and sensitivity of biomass to a 10% mortality increase, and the
# consequences of species removals. Step sizes and the number of extinction sequences follow
# what Dataset S1 shows was done (lab/wiki/model/simulator-reconstruction.md).

function Base.copy(m::FoodWebModel)
    return FoodWebModel(m.web, m.params, m.exponents, m.trophic_level, m.mass, m.basal, m.growth, m.metabolism,
                        m.base_mortality, copy(m.mortality), m.efficiency, m.attack, m.handling, m.resources,
                        copy(m.alive), copy(m.w), m.threshold, zeros(length(m.Bp)), zeros(length(m.Bp)), zeros(length(m.Bp)),
                        m.link_attack, m.link_ah)
end

"""
    PerturbationOptions(; press = 0.1, global_start = 1.001, step = 0.1, coarse = 10,
                          max_steps = 2000, sequences = 10,
                          steady = SteadyStateOptions(tmax = 1e4, maxiters = 2 * 10^5, stall = 1e3))

Settings of the press experiments: the size of the mortality press for resistance (`press`,
10%), the first global mortality multiplier and the grid step of the tolerance searches
(`global_start`, `step`), the first jump of the tolerance search in grid steps (`coarse`; it
doubles while no species dies and halves when one does, down to single grid steps), the largest number of grid steps tried, the number of random extinction sequences
for tolerance to extinctions (`sequences`; the SI says 100, Dataset S1 shows 10), and the
steady-state options used after each perturbation. After a perturbation the integration is
limited to `tmax = 10⁴` (assembly uses T_SAFE = 10⁶), and stops once no species has died for
1000 time units: such a community is oscillating, and its state at that time is used. SI §10 does not say how
non-stationary perturbed states were treated.
"""
Base.@kwdef struct PerturbationOptions
    press::Float64 = 0.1
    global_start::Float64 = 1.001
    step::Float64 = 0.1
    coarse::Int = 10
    max_steps::Int = 2000
    sequences::Int = 10
    steady::SteadyStateOptions = SteadyStateOptions(tmax = 1e4, maxiters = 2 * 10^5, stall = 1e3)
end

# Move `model` (already modified) to its new steady state from `B`. A Newton step from the old
# equilibrium is accepted when it lands on a stable point with all survivors above threshold
# (`newton = true`); otherwise the model is integrated, which also handles extinctions.
function equilibrate!(model::FoodWebModel, B::Vector{Float64}, steady::SteadyStateOptions; newton::Bool = true)
    if newton
        Bn = copy(B)
        alive0 = copy(model.alive)
        if polish!(model, Bn)
            copyto!(B, Bn)
            return true
        end
        copyto!(model.alive, alive0)
    end
    Bs, ok = settle!(model, B; options = steady)
    copyto!(B, Bs)
    return ok
end

# Smallest grid index k ≥ k0 at which raising mortality by `factor(k)` (applied by `apply!`)
# drives some species extinct. Each level is reached from the last extinction-free state (by
# Newton when the equilibrium stays stable, otherwise by integration). Jumps start at
# `opts.coarse` grid steps and double while no species dies; when a jump kills one it is
# halved and retried, so the deciding step is always a single grid increment, as in the
# paper's step-by-step increase. This assumes a jump that causes no extinction would cause
# none if taken in smaller steps.
function first_extinction(c::Community, apply!::Function, k0::Integer, opts::PerturbationOptions)
    model = copy(c.model)
    B = copy(c.B)
    n0 = richness(model)
    function attempt(k)
        trial = copy(model)
        Bt = copy(B)
        apply!(trial, k)
        equilibrate!(trial, Bt, opts.steady)
        return richness(trial) == n0, trial, Bt
    end
    ok, trial, Bt = attempt(k0)
    ok || return k0
    model, B = trial, Bt
    lo, jump = k0, max(opts.coarse, 1)
    while lo < opts.max_steps
        hi = min(lo + jump, opts.max_steps)
        ok, trial, Bt = attempt(hi)
        if ok
            model, B, lo = trial, Bt, hi
            jump *= 2
        elseif hi == lo + 1
            return hi
        else
            jump = max((hi - lo) ÷ 2, 1)
        end
    end
    return nothing
end

"""
    global_mortality_tolerance(c::Community; options = PerturbationOptions())

The raw column `TM_p`: the mortality multiplier `global_start + step k` (k = 0, 1, …) applied
to every species at which the first extinction occurs, raising it step by step from the
equilibrium. Dataset S1 records the multiplier itself (1.001, 1.101, …); SI eq. 17 writes it
as a relative increase. `NaN` if no extinction within `max_steps`.
"""
function global_mortality_tolerance(c::Community; options::PerturbationOptions = PerturbationOptions())
    f(k) = options.global_start + options.step * k
    k = first_extinction(c, (m, k) -> set_mortality!(m, f(k)), 0, options)
    return k === nothing ? NaN : f(k)
end

"""
    local_mortality_tolerances(c::Community; options = PerturbationOptions())

For each living species in turn, the smallest relative increase `step k` (k = 1, 2, …) in its
own mortality that causes any extinction. Returns a vector in the order of
`findall(c.model.alive)`; the raw columns `minTM_i` and `av_TM_i` are its minimum and mean.
"""
function local_mortality_tolerances(c::Community; options::PerturbationOptions = PerturbationOptions())
    out = Float64[]
    for i in findall(c.model.alive)
        k = first_extinction(c, (m, k) -> (m.mortality[i] = 1 + options.step * k), 1, options)
        push!(out, k === nothing ? NaN : options.step * k)
    end
    return out
end

"""
    press_response(c::Community, apply!; options = PerturbationOptions())

Apply `apply!(model)` to a copy of the community, move it to its new steady state, and return
`(relative_change, l1_change, extinctions)`: `(B₀ − B′)/B₀` for total biomass, the L1 norm of
the species-wise change divided by B₀, and the number of species lost.
"""
function press_response(c::Community, apply!::Function; options::PerturbationOptions = PerturbationOptions(),
                        newton::Bool = true)
    model = copy(c.model)
    B = copy(c.B)
    n0 = richness(model)
    apply!(model, B)
    equilibrate!(model, B, options.steady; newton = newton)
    B0 = sum(c.B)
    return ((B0 - sum(B)) / B0, sum(abs, c.B .- B) / B0, n0 - richness(model))
end

"""
    removal_responses(c::Community; options = PerturbationOptions())

Remove each living species in turn from the equilibrium community and integrate to the new
steady state. Returns a vector of `(survivors, relative_change, l1_change)`, where `survivors`
counts the species alive afterwards (at most N − 1).
"""
function removal_responses(c::Community; options::PerturbationOptions = PerturbationOptions())
    out = NamedTuple{(:survivors, :relative, :l1),Tuple{Int,Float64,Float64}}[]
    B0 = sum(c.B)
    for i in findall(c.model.alive)
        model = copy(c.model)
        B = copy(c.B)
        kill!(model, B, i)
        equilibrate!(model, B, options.steady; newton = false)
        push!(out, (survivors = richness(model), relative = (B0 - sum(B)) / B0, l1 = sum(abs, c.B .- B) / B0))
    end
    return out
end

"""
    extinction_tolerance([rng], c::Community; options = PerturbationOptions())

Robustness (Dunne et al. 2002): the number of random primary removals, each followed by
integration to steady state, needed before at most half of the original species remain,
averaged over `options.sequences` random sequences. The raw column `Ext_half`.
"""
function extinction_tolerance(rng::AbstractRNG, c::Community; options::PerturbationOptions = PerturbationOptions())
    n0 = richness(c)
    total = 0
    for _ in 1:options.sequences
        model = copy(c.model)
        B = copy(c.B)
        removed = 0
        while richness(model) > n0 / 2
            i = rand(rng, findall(model.alive))
            kill!(model, B, i)
            equilibrate!(model, B, options.steady; newton = false)
            removed += 1
        end
        total += removed
    end
    return total / options.sequences
end
