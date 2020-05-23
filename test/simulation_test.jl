using BeastUtils.Simulation, BeastUtils.TreeUtils
using LinearAlgebra, Test, PhyloNetworks

import Random
seed = 666
Random.seed!(seed)


function check_simulation(tsm::Simulation.TraitSimulationModel,
                            μ_true::Vector{Float64},
                            V_true::Matrix{Float64};
                            reps::Int = 200_000,
                            μ_tol::Float64 = 1e-1,
                            v_tol::Float64 = 1e-1)

    V = zeros(n * p, n * p)
    μ = zeros(n * p)


    for i = 1:reps
        y = vec(Simulation.simulate(tsm))
        μ .+= y
        V .+= y * y'
    end

    μ ./= reps

    V ./= reps
    V = V - μ * μ'

    V_diff = abs.(V - V_true)
    μ_diff = abs.(μ - μ_true)

    v_max = maximum(V_diff)
    μ_max = maximum(μ_diff)



    @test v_max < v_tol
    @test μ_max < μ_tol
end

n = 5
p = 3
taxa = ["taxon$i" for i = 1:n]



## Tree only



tree = TreeUtils.rtree(n, labels = taxa, keep_order = false, ultrametric = false)
Σ = randn(p, p)
Σ = Σ * Σ' # needs to be positive definite
tdm = Simulation.TreeDiffusionModel(tree, Σ)
tsm = Simulation.TraitSimulationModel(taxa, tdm)

μ_true = zeros(n * p)
Ψ = TreeUtils.vcv(tree, taxa)
V_true = kron(Σ, Ψ)

check_simulation(tsm, μ_true, V_true)


## Residual variance

tree = TreeUtils.rtree(n, labels = taxa, keep_order = false, ultrametric = false)


Σ = randn(p, p)
Γ = randn(p, p)


Σ = Σ * Σ' # needs to be positive definite
Γ = Γ * Γ' # needs to be positive definite


tdm = Simulation.TreeDiffusionModel(tree, Σ)
rvm = Simulation.ResidualVarianceModel(Γ)
tsm = Simulation.TraitSimulationModel(taxa, tdm, rvm)


μ_true = zeros(n * p)
Ψ = TreeUtils.vcv(tdm.tree, taxa)
V_true = kron(Σ, Ψ) + kron(Γ, Diagonal(ones(n)))

check_simulation(tsm, μ_true, V_true)


## Latent factor model

k = 2


tree = TreeUtils.rtree(n, labels = taxa, keep_order = false, ultrametric = false)


Σ = Diagonal(ones(k))
Λ = Diagonal(rand(p))
L = randn(k, p)




tdm = Simulation.TreeDiffusionModel(tree, Σ)
lfm = Simulation.LatentFactorModel(L, Λ)
tsm = Simulation.TraitSimulationModel(taxa, tdm, lfm)


μ_true = zeros(n * p)
Ψ = TreeUtils.vcv(tdm.tree, taxa)
V_true = kron(L' * L, Ψ) + kron(Λ, Diagonal(ones(n)))

check_simulation(tsm, μ_true, V_true)
