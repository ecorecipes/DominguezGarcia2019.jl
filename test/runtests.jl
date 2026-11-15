using Test
using EcologicalResilience
using CSV, DataFrames
using LinearAlgebra: Diagonal, I, eigvals, norm
using Random: Xoshiro

const REF = joinpath(@__DIR__, "reference")

@testset "EcologicalResilience" begin
    include("test_niche.jl")
    include("test_bioenergetic.jl")
    include("test_local_metrics.jl")
    include("test_assembly.jl")
    include("test_analysis.jl")
    include("test_matrices.jl")
    include("test_timeseries.jl")
    include("test_docstrings.jl")
end
