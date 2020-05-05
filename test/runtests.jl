using Test, SafeTestsets


@time @safetestset "Data simulation test" begin include("diffusionSimulation_test.jl") end
@time @safetestset "RTrees Test" begin include("RTrees_test.jl") end
