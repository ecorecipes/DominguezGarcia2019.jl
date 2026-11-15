# Readers for published interaction datasets (local copies; see lab/data/<dataset>/README.md for
# provenance). Each returns typed matrices whose metadata records the conversion applied.

"""
    read_kinlock(path; control = :monoculture, complete = true)

Plant interaction networks of Kinlock (2019, Dryad doi:10.5061/dryad.1sm06sp, CC0): for each
study, a matrix of relative interaction indices (RII ∈ [−1, 1]; entry (i, j) the effect of
species j on species i, negative for competition), against monoculture controls
(`control = :monoculture`, column `NetworkMonoCtrl`, diagonal 0 by construction) or plants
grown alone (`:alone`, `NetworkTrueCtrl`). With `complete = true` only networks without
missing entries are returned. Returns a vector of `(study, species, rii, habitat)`.
"""
function read_kinlock(path::AbstractString; control::Symbol = :monoculture, complete::Bool = true)
    col = control === :monoculture ? "NetworkMonoCtrl" : control === :alone ? "NetworkTrueCtrl" :
          throw(ArgumentError("control must be :monoculture or :alone, got $control"))
    df = CSV.read(path, DataFrame; stringtype = String)
    out = NamedTuple{(:study, :species, :rii, :habitat),Tuple{String,Vector{String},Matrix{Float64},String}}[]
    for r in eachrow(df)
        s = r[col]
        (s === missing || s == "NA") && continue
        vals = [v == "NA" ? NaN : parse(Float64, v) for v in strip.(split(s, ","))]
        n = Int(r.SpeciesRichness)
        length(vals) == n^2 || continue
        M = permutedims(reshape(vals, n, n))       # stored row by row
        complete && any(isnan, M) && continue
        sp = String.(strip.(split(r.Species, ",")))
        push!(out, (study = r.Study, species = sp, rii = M, habitat = string(r.Habitat)))
    end
    return out
end

"""
    alpha_from_rii(net; scale = 1.0)

An [`AlphaMatrix`](@ref) from a Kinlock network: off-diagonal per-capita effects
`scale × RII` relative to self-regulation −1. RII values are relative changes in performance,
not per-capita coefficients, so this is an assumption, recorded in the metadata; the growth
rates are unknown and set to 1, which affects equilibria but not the feasibility domain.
"""
function alpha_from_rii(net; scale::Real = 1.0)
    A = scale .* net.rii
    A[diagind(A)] .= -1.0
    meta = MatrixMetadata(; taxa = net.species, state_variables = "plant performance (study-specific)",
                          units = "relative to self-regulation", spatial_bounds = net.habitat,
                          derivation = "Kinlock 2019 RII as per-capita effects × $scale, diagonal −1; r unknown, set to 1",
                          licence = "CC0-1.0")
    return AlphaMatrix(A, ones(length(net.species)); meta = meta)
end

"""
    read_adler(path; normalise = true)

Complete competition-coefficient matrices from the database of Adler et al. (2018, Ecology
Letters, Dryad doi:10.5061/dryad.q5mg97b, CC0): rows grouped by study, response variable and
treatment, kept when every focal × competitor pair (including intraspecific) is present.
Signs are put in the Lotka–Volterra convention (negative = competition) using the column
"Negative coefs mean competition" (`Yes`/`No`; other entries skip the group). With
`normalise = true` each row is divided by the absolute intraspecific coefficient, so the
diagonal is −1. Growth rates are not in the database (set to 1). Returns a vector of
`(study, response, treatment, alpha::AlphaMatrix)`.
"""
function read_adler(path::AbstractString; normalise::Bool = true)
    df = CSV.read(path, DataFrame; stringtype = String)
    rename!(df, Dict("Paper Key #" => "key", "First Author" => "author", "Response variable" => "response",
                     "competition coefficient" => "coef", "Negative coefs mean competition" => "negcomp"))
    df.treatment = coalesce.(df.treatment, "")
    out = NamedTuple{(:study, :response, :treatment, :alpha),Tuple{String,String,String,AlphaMatrix}}[]
    for g in groupby(df, [:key, :author, :Year, :response, :treatment])
        nc = unique(skipmissing(g.negcomp))
        length(nc) == 1 && nc[1] in ("Yes", "No") || continue
        sp = sort(unique(vcat(g[!, "Focal.Species"], g[!, "Comp.Species"])))
        n = length(sp)
        n ≥ 2 || continue
        pos = Dict(s => k for (k, s) in enumerate(sp))
        A = fill(NaN, n, n)
        for r in eachrow(g)
            ismissing(r.coef) && continue
            c = r.coef isa Number ? Float64(r.coef) : tryparse(Float64, string(r.coef))
            c === nothing && continue
            A[pos[r["Focal.Species"]], pos[r["Comp.Species"]]] = nc[1] == "Yes" ? c : -c
        end
        any(isnan, A) && continue
        if normalise
            d = abs.(diag(A))
            any(iszero, d) && continue
            A = A ./ d
        end
        note = string("Adler et al. 2018 coefficients, sign to gLV convention", normalise ? ", rows / |αᵢᵢ|" : "",
                      "; r unknown, set to 1")
        meta = MatrixMetadata(; taxa = sp, state_variables = string(g.response[1]), derivation = note,
                              licence = "CC0-1.0")
        study = string(g.author[1], " ", g.Year[1])
        push!(out, (study = study, response = string(g.response[1]), treatment = string(g.treatment[1]),
                    alpha = AlphaMatrix(A, ones(n); meta = meta)))
    end
    return out
end

"""
    read_stein(dir)

The gut-microbiota generalised Lotka–Volterra model and data of Stein et al. (2013, PLoS
Comput Biol 9:e1003388): `dir` holds `derived/stein2013_glv.csv` (the published fit, sheet
MmuE of their Dataset S1) and `clv_data_stein/raw_data.csv` (the time series, MIT-licensed
copy from the cLV repository). Returns `(alpha, susceptibility, data)`: the published
[`AlphaMatrix`](@ref) (rows and columns in the sheet's order), the clindamycin susceptibilities
ε (a pulse `u(t)` adds `εᵢ u(t)` to the per-capita growth of species i), and a long `DataFrame`
with columns `mouse`, `group`, `time`, `genus`, `abundance` (10¹¹ 16S copies per cm³) and
`clindamycin`.
"""
function read_stein(dir::AbstractString)
    fit = CSV.read(joinpath(dir, "derived", "stein2013_glv.csv"), DataFrame; stringtype = String)
    genera = fit.genus
    A = Matrix{Float64}(fit[:, genera])
    meta = MatrixMetadata(; taxa = genera, state_variables = "16S rRNA gene copies (10¹¹ per cm³ of caecal content)",
                          units = "per day", spatial_bounds = "gnotobiotic mice",
                          derivation = "Stein et al. 2013, Dataset S1 sheet MmuE (published regularised fit)",
                          licence = "CC BY 4.0 (article)")
    alpha = AlphaMatrix(A, fit.growth; meta = meta)
    raw = CSV.read(joinpath(dir, "clv_data_stein", "raw_data.csv"), DataFrame; header = false, stringtype = String)
    lab = String.(raw[:, 1])
    row(name) = Vector(raw[findfirst(==(name), lab), 2:end])
    num(v) = [x isa Number ? Float64(x) : parse(Float64, string(x)) for x in v]
    group, mouse, time, clin = Int.(num(row("Population"))), Int.(num(row("ID"))), num(row("time (in d)")),
                               num(row("Clindamycin signal"))
    data = DataFrame(mouse = Int[], group = Int[], time = Float64[], genus = String[], abundance = Float64[],
                     clindamycin = Float64[])
    for g in genera
        v = num(row(g))
        append!(data, DataFrame(mouse = mouse, group = group, time = time, genus = g, abundance = v, clindamycin = clin))
    end
    return (alpha = alpha, susceptibility = Vector{Float64}(fit.susceptibility), data = data)
end
