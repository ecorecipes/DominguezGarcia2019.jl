using Documenter
using DominguezGarcia2019

tutorials = sort(filter(f -> endswith(f, ".md"), readdir(joinpath(@__DIR__, "src", "tutorials"))))
titles = Dict("01_reproducing_the_analysis.md" => "Reproducing the analysis",
              "02_one_community.md" => "One food web and its 27 metrics",
              "03_simulated_ensemble.md" => "End to end: a simulated ensemble",
              "04_early_response.md" => "The early response to pulse: an indexing artefact",
              "05_observed_food_webs.md" => "Observed food webs (GATEWAy)",
              "06_robustness.md" => "Robustness to model settings")

# Source and edit links need a git remote; outside a git checkout, build without them.
ingit = success(pipeline(`git -C $(@__DIR__) rev-parse`; stdout = devnull, stderr = devnull))

makedocs(;
         modules = [DominguezGarcia2019],
         authors = "Simon Frost",
         sitename = "DominguezGarcia2019.jl",
         checkdocs = :exports,
         remotes = ingit ? Dict() : nothing,
         format = Documenter.HTML(;
                                  edit_link = ingit ? "main" : nothing,
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://ecorecipes.github.io/DominguezGarcia2019.jl",
                                  size_threshold = 400 * 1024),
         pages = ["Home" => "index.md",
                  "Tutorials" => [get(titles, f, f) => "tutorials/$f" for f in tutorials],
                  "API Reference" => ["Datasets and conventions" => "api/metrics.md",
                                      "Analysis" => "api/analysis.md"]])

if get(ENV, "CI", "false") == "true"
    deploydocs(; repo = "github.com/ecorecipes/DominguezGarcia2019.jl.git", push_preview = true)
end
