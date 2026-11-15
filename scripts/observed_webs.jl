# Stability metrics for observed food webs (GATEWAy v1.0) under the bioenergetic model: the
# observed topology replaces the niche model; each web is run REPLICATES times with its own draw
# of allometric exponents and initial biomasses.
#
#   julia -t auto --gcthreads=1 --project scripts/observed_webs.jl GATEWAY_CSV OUT.csv [REPLICATES MAXNODES SEED MASSES NODES]
#
# NODES = trophic (default) aggregates taxa with identical prey and predators into trophic species
# (`trophic_species`), as Williams & Martinez (2000) did when testing the niche model; NODES = taxa
# keeps the published nodes. Unaggregated, identical taxa get identical parameters and survive
# together, so the equilibrium is neutral (Jacobian eigenvalue 0): in a run on taxa, over a quarter
# of the converged communities had Rinf = 0. MAXNODES applies after aggregation.
#
# MASSES = trophic (default) uses the model's masses Z^(L−1) on the observed topology; MASSES =
# measured uses GATEWAy body masses relative to the producers. With measured masses most
# consumers starve: in this parameterisation maximum intake per unit biomass scales as m^−0.4
# and metabolism as m^−0.25, so realistic mass ratios (median 3×10⁵) are not viable
# (lab/wiki/analysis/observed-webs.md). A web is used if its trophic levels are defined, it has
# at most MAXNODES nodes, and (for measured masses) every consumer has a mass. Non-living
# resources (detritus) are basal nodes with logistic renewal, an approximation recorded in the
# output (`nonliving`). There is no survival filter: the observed
# web is not a random draw to be accepted or rejected. Communities that do not reach a stable
# equilibrium, or keep fewer than 3 species, are recorded with `status` and no metrics. `merged`
# counts neutral top consumers merged (see `merge_neutral!`). Rows are
# appended as they finish, so an interrupted run resumes.
using EcologicalResilience
using CSV, DataFrames
using LinearAlgebra: BLAS
using Random: Xoshiro

BLAS.set_num_threads(1)

src, out = ARGS[1], ARGS[2]
reps = length(ARGS) ≥ 3 ? parse(Int, ARGS[3]) : 5
maxnodes = length(ARGS) ≥ 4 ? parse(Int, ARGS[4]) : 120
seed = length(ARGS) ≥ 5 ? parse(Int, ARGS[5]) : 20190000
masses = length(ARGS) ≥ 6 ? ARGS[6] : "trophic"
masses in ("trophic", "measured") || error("MASSES must be trophic or measured")
nodes = length(ARGS) ≥ 7 ? ARGS[7] : "trophic"
nodes in ("trophic", "taxa") || error("NODES must be trophic or taxa")

taxa_webs = read_gateway(src)
webs = nodes == "trophic" ? trophic_species.(taxa_webs) : taxa_webs
usable = [k for (k, w) in enumerate(webs)
          if richness(w.web) ≤ maxnodes && trophic_levels(w.web) !== nothing &&
             (masses == "trophic" || relative_masses(w) !== nothing)]
jobs = [(k, r) for k in usable for r in 1:reps]
id(k, r) = "GW$(k)_$(r)"
done = isfile(out) ? Set(CSV.read(out, DataFrame; select = [:network]).network) : Set{String}()
todo = [j for j in jobs if id(j...) ∉ done]
println("$(length(usable)) usable webs, $(length(jobs)) jobs, $(length(todo)) to do, $(Threads.nthreads()) threads")
lk = ReentrantLock()
header = !isfile(out)
params = BioenergeticParameters()
empty_metrics = NamedTuple{Tuple(RAW_COLUMNS)}(ntuple(_ -> NaN, length(RAW_COLUMNS)))

# Surviving consumers with the same surviving prey, the same mass and no surviving predators are
# dynamically identical and exactly neutral (their growth is linear in their own biomass), even
# after aggregation into trophic species: extinctions can remove the links that told them apart.
# Merging them into one species with their summed biomass leaves the dynamics unchanged. Identical
# species that do have predators are not neutral: predation with Hill exponent 1 + q favours
# the rarer one.
function merge_neutral!(model, B)
    A, alive = model.web.A, model.alive
    top = [i for i in findall(alive) if !model.basal[i] && !any(alive[k] && A[k, i] for k in eachindex(alive))]
    merged = 0
    for (a, i) in enumerate(top), j in top[(a + 1):end]
        (alive[i] && alive[j] && model.mass[i] == model.mass[j]) || continue
        if all(A[i, k] == A[j, k] for k in eachindex(alive) if alive[k])
            B[i] += B[j]
            kill!(model, B, j)
            merged += 1
        end
    end
    return merged
end

Threads.@threads :dynamic for (k, r) in todo
    w = webs[k]
    rng = Xoshiro(seed + 1000k + r)
    S = richness(w.web)
    model = FoodWebModel(w.web, params, draw_exponents(rng, params);
                         mass = masses == "measured" ? relative_masses(w) : nothing)
    B, ok = settle!(model, 0.05 .+ 0.95 .* rand(rng, S); options = SteadyStateOptions(tmax = 1e5))
    nd = richness(model)
    iso = isolated_plants(model)
    isempty(iso) || kill!(model, B, iso)
    merged = ok ? merge_neutral!(model, B) : 0
    status = !ok ? "not_converged" : richness(model) < 3 ? "collapsed" : "ok"
    rec = empty_metrics
    t = @elapsed if status == "ok"
        try
            # :si: the published median weights assume niche-ordered species, which an observed web lacks
            rec = community_record(rng, Community(model, B, S, nd); convention = :si)
        catch err
            status = "error: " * first(sprint(showerror, err), 80)
        end
        # a zero eigenvalue left by dynamically identical survivors: no meaningful metrics
        if status == "ok" && !(rec.Rinf > 1e-9)
            status, rec = "neutral", empty_metrics
        end
    end
    info = (network = id(k, r), web = w.name, ecosystem = w.ecosystem, replicate = r, masses = masses, node_type = nodes,
            taxa = richness(taxa_webs[k].web), nodes = S, nonliving = count(.!w.living), survivors_dynamic = nd,
            status = status, merged = merged)
    row = DataFrame([merge(info, rec, (seconds = round(t; digits = 2),))])
    lock(lk) do
        CSV.write(out, row; append = !header)
        global header = false
        println("$(id(k, r)) $(w.ecosystem) S=$S survivors=$(richness(model)) $status $(round(t; digits = 1)) s")
        flush(stdout)
    end
end
