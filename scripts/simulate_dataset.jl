# Simulate a dataset of stability metrics in the layout of the raw Dataset S1 of
# Domínguez-García et al. (2019): one row per assembled community, all raw columns.
#
#   julia -t auto --gcthreads=1 --project scripts/simulate_dataset.jl OUT.csv NCOMMUNITIES [MINS MAXS SEED [KEY=VALUE ...]]
#
# KEY=VALUE pairs after the positional arguments override fields of BioenergeticParameters,
# e.g. `q=0.6 h0=0.2 plant_metabolism=true`, for sensitivity analyses.
#
# --gcthreads=1: Julia 1.12's parallel GC segfaulted (gc_sweep_pool_parallel) after about 800
# communities; a run with --check-bounds=yes found no out-of-bounds access in this package.
#
# Community k uses the random stream Xoshiro(SEED + k), draws its initial richness uniformly
# from MINS:MAXS (SI: 5–115; Dataset S1: 5–110) and redraws webs until one is viable, as SI §9
# describes. Rows are appended as they finish, so an interrupted run resumes where it stopped.
# The raw columns follow the published convention (`community_record`); `Amax_si`, `tmax_si`,
# `MRinf_si`, `MR0_si`, `MAmax_si` and `Mtmax_si` give the same six metrics as the SI defines them
# (`amplification_metrics(c; convention = :si)`), and `RM_p_si` and `L1DBRM_p_si` the global press
# as the SI defines it (every mortality raised by 10%, `global_press_response`).
using DominguezGarcia2019
using FoodWebDynamics
using CSV, DataFrames
using LinearAlgebra: BLAS
using Random: Xoshiro

BLAS.set_num_threads(1)

out = ARGS[1]
n = parse(Int, ARGS[2])
mins = length(ARGS) ≥ 3 ? parse(Int, ARGS[3]) : 5
maxs = length(ARGS) ≥ 4 ? parse(Int, ARGS[4]) : 110
seed = length(ARGS) ≥ 5 ? parse(Int, ARGS[5]) : 20191217
overrides = map(ARGS[6:end]) do a
    k, v = Symbol.(split(a, "="; limit = 2))
    T = fieldtype(BioenergeticParameters, k)
    k => parse(T, String(v))
end
params = BioenergeticParameters(; overrides...)
isempty(overrides) || println("model settings: ", join(("$k = $v" for (k, v) in overrides), ", "))

done = isfile(out) ? Set(CSV.read(out, DataFrame; select = [:network]).network) : Set{String}()
todo = [k for k in 1:n if "Sim$(seed + k)" ∉ done]
println("$(length(done)) done, $(length(todo)) to do, $(Threads.nthreads()) threads")
lock_ = ReentrantLock()
header = !isfile(out)

Threads.@threads :dynamic for k in todo
    rng = Xoshiro(seed + k)
    S = rand(rng, mins:maxs)
    attempts = Dict{Symbol,Int}()
    c = nothing
    while c === nothing
        c, status = assemble_community(rng, S; params = params)
        attempts[status] = get(attempts, status, 0) + 1
    end
    rec = nothing
    t = @elapsed try
        si = amplification_metrics(c; convention = :si)
        rel, l1, _ = DominguezGarcia2019.global_press_response(c)
        rec = merge(community_record(rng, c),
                    NamedTuple{Symbol.(string.(keys(si), "_si"))}(values(si)),
                    (RM_p_si = rel, L1DBRM_p_si = l1))
    catch err                       # log and skip; a rerun retries it (and fails the same way)
        lock(lock_) do
            println("Sim$(seed + k): S=$S ERROR ", first(sprint(showerror, err), 200))
            flush(stdout)
        end
    end
    rec === nothing && continue
    row = DataFrame([merge((network = "Sim$(seed + k)",), rec,
                           (attempts = sum(values(attempts)), seconds = round(t; digits = 2)))])
    lock(lock_) do
        CSV.write(out, row; append = !header)
        global header = false
        println("Sim$(seed + k): S=$S N=$(rec.N) $(round(t; digits = 1)) s")
        flush(stdout)
    end
end
