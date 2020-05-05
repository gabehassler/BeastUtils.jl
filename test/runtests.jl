using Test, SafeTestsets

@time @safetestset "Data simulation test" begin include("diffusionSimulation_test.jl") end
