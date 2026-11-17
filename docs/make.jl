using Documenter
using EcologicalResilience

tutorials = sort(filter(f -> endswith(f, ".md"), readdir(joinpath(@__DIR__, "src", "tutorials"))))
titles = Dict("01_reproducing_the_analysis.md" => "Reproducing the analysis",
              "02_one_community.md" => "One food web and its 27 metrics",
              "03_simulated_ensemble.md" => "End to end: a simulated ensemble",
              "04_early_response.md" => "The early response to pulse: an indexing artefact")

# Source and edit links need a git remote; outside a git checkout, build without them.
ingit = success(pipeline(`git -C $(@__DIR__) rev-parse`; stdout = devnull, stderr = devnull))

makedocs(;
         modules = [EcologicalResilience],
         authors = "Simon Frost",
         sitename = "EcologicalResilience.jl",
         checkdocs = :exports,
         remotes = ingit ? Dict() : nothing,
         format = Documenter.HTML(;
                                  edit_link = ingit ? "main" : nothing,
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://ecorecipes.github.io/EcologicalResilience.jl",
                                  size_threshold = 400 * 1024),
         pages = ["Home" => "index.md",
                  "Tutorials" => [get(titles, f, f) => "tutorials/$f" for f in tutorials],
                  "API Reference" => ["Stability metrics" => "api/metrics.md",
                                      "Analysis" => "api/analysis.md",
                                      "Interaction matrices" => "api/matrices.md"]])

if get(ENV, "CI", "false") == "true"
    deploydocs(; repo = "github.com/ecorecipes/EcologicalResilience.jl.git", push_preview = true)
end
