# Recompute the global-press columns (RM_p, L1DBRM_p) of an ensemble written by
# simulate_dataset.jl with the published press (`published_press_factors`), without rerunning the
# other experiments.
#
#   julia -t auto --gcthreads=1 --project scripts/recompute_global_press.jl ENSEMBLE.csv OUT.csv [MINS MAXS]
#
# Community `Sim<s>` is re-assembled from Xoshiro(s) as simulate_dataset.jl did (default model
# settings; MINS:MAXS as in that run, 5:110 by default). Its resilience is checked against the
# stored Rinf, so a community that does not reproduce is flagged rather than silently replaced.
# OUT.csv gets `network`, `RM_p`, `L1DBRM_p`, `Rinf_match`; rows are appended as they finish, so an
# interrupted run resumes.
using DominguezGarcia2019
using FoodWebDynamics
using CSV, DataFrames
using LinearAlgebra: BLAS
using Random: Xoshiro

BLAS.set_num_threads(1)

ens = CSV.read(ARGS[1], DataFrame; select = [:network, :Rinf])
out = ARGS[2]
mins = length(ARGS) ≥ 3 ? parse(Int, ARGS[3]) : 5
maxs = length(ARGS) ≥ 4 ? parse(Int, ARGS[4]) : 110

done = isfile(out) ? Set(CSV.read(out, DataFrame; select = [:network]).network) : Set{String}()
todo = [r for r in eachrow(ens) if r.network ∉ done]
println("$(length(done)) done, $(length(todo)) to do, $(Threads.nthreads()) threads")
lock_ = ReentrantLock()
header = !isfile(out)

Threads.@threads :dynamic for r in todo
    rng = Xoshiro(parse(Int, r.network[4:end]))
    S = rand(rng, mins:maxs)
    c = nothing
    while c === nothing
        c, _ = assemble_community(rng, S)
    end
    row = try
        rel, l1, _ = press_response(c, published_press_factors(c))
        (network = r.network, RM_p = rel, L1DBRM_p = l1,
         Rinf_match = isapprox(resilience(community_jacobian(c)), r.Rinf; rtol = 1e-6))
    catch err
        lock(lock_) do
            println("$(r.network): ERROR ", first(sprint(showerror, err), 200))
            flush(stdout)
        end
        nothing
    end
    row === nothing && continue
    lock(lock_) do
        CSV.write(out, DataFrame([row]); append = !header)
        global header = false
        row.Rinf_match || (println("$(r.network): Rinf does not match the stored row"); flush(stdout))
    end
end
println("finished")
