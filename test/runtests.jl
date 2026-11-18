using Test
using DominguezGarcia2019
using FoodWebDynamics
using CSV, DataFrames
using LinearAlgebra: Diagonal, I, eigvals, norm
using Random: Xoshiro

const REF = joinpath(@__DIR__, "reference")

@testset "DominguezGarcia2019" begin
    include("test_conventions.jl")
    include("test_assembly.jl")
    include("test_analysis.jl")
    include("test_docstrings.jl")
end
