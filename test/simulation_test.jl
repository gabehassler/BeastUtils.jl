using BeastUtils.Simulation, BeastUtils.TreeUtils, BeastUtils.SlowLikelihoods,
      BeastUtils.MatrixUtils
using LinearAlgebra, Test, PhyloNetworks

import Random
seed =666
Random.seed!(seed)

const var_diff_tol = 1e-10


function check_simulation(tsm::Simulation.TraitSimulationModel,
                            μ_true::Vector{Float64},
                            V_true::Matrix{Float64};
                            reps::Int = 100_000,
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

    @show v_max
    @show μ_max

    return v_max < v_tol && μ_max < μ_tol
end

n = 5
p = 3
taxa = ["taxon$i" for i = 1:n]



## Tree only

# Known tree
newick = "((A:1.0,(B:0.5,C:0.8):0.1):2.0,(D:0.5,E:2.0):1.5);"
tree = readTopology(newick)
Σ = randpd(p, rescale=true)
Ψ1 = [3.0 2.0 2.0 0.0 0.0;
      2.0 2.6 2.1 0.0 0.0;
      2.0 2.1 2.9 0.0 0.0;
      0.0 0.0 0.0 2.0 1.5;
      0.0 0.0 0.0 1.5 3.5]

taxa1 = ["A", "B", "C", "D", "E"]
μ_true = zeros(n * p)
V_true = kron(Σ, Ψ1)
tdm = Simulation.TreeDiffusionModel(tree, Σ)
tsm = Simulation.TraitSimulationModel(taxa1, tdm)

model_params = SlowLikelihoods.sim_to_params(tsm, pss=Inf)
V_true2 = var(model_params, tree, taxa1)

@test maximum(abs.(V_true - V_true2)) < var_diff_tol


@test check_simulation(tsm, μ_true, V_true)

perm = [3, 2, 5, 4, 1] # scramble order
taxa2 = taxa1[perm]
Ψ2 = Ψ1[perm, perm]

μ_true = zeros(n * p)
V_true = kron(Σ, Ψ2)
tdm = Simulation.TreeDiffusionModel(tree, Σ)
tsm = Simulation.TraitSimulationModel(taxa2, tdm)

model_params = SlowLikelihoods.sim_to_params(tsm, pss=Inf)
V_true2 = var(model_params, tree, taxa2)
@test maximum(abs.(V_true - V_true2)) < var_diff_tol

@test check_simulation(tsm, μ_true, V_true)

# Random Tree
tree = TreeUtils.rtree(n, labels = taxa, keep_order = false, ultrametric = true)
standardize_height!(tree)
Σ = randpd(p, rescale=true)
tdm = Simulation.TreeDiffusionModel(tree, Σ)
tsm = Simulation.TraitSimulationModel(taxa, tdm)

μ_true = zeros(n * p)
Ψ = TreeUtils.vcv(tree, taxa)
V_true = kron(Σ, Ψ)

model_params = SlowLikelihoods.sim_to_params(tsm, pss=Inf)
V_true2 = var(model_params, tree, taxa)
@test maximum(abs.(V_true - V_true2)) < var_diff_tol

@test check_simulation(tsm, μ_true, V_true)

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

model_params = SlowLikelihoods.sim_to_params(tsm, pss=Inf)
V_true2 = var(model_params, tree, taxa)
@test maximum(abs.(V_true - V_true2)) < var_diff_tol

@test check_simulation(tsm, μ_true, V_true)


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

@test check_simulation(tsm, μ_true, V_true)
