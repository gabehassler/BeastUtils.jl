using BeastUtils.DiffusionSimulation, BeastUtils.RTrees
using LinearAlgebra

import Random
seed = 666
Random.seed!(seed)

n = 5
p = 4

Σ = randn(p, p)
Γ = randn(p, p)

Σ = Diagonal(ones(p))
Γ = Diagonal(ones(p))

Σ = Σ * Σ' # needs to be positive definite
Γ = Γ * Γ' # needs to be positive definite

taxa = ["taxon$i" for i = 1:n]
tree = RTrees.rtree(taxa, seed)
newick = RTrees.make_newick(tree)

tdm = DiffusionSimulation.TreeDiffusionModel(newick, Σ, Γ)

V = zeros(n * p, n * p)
μ = zeros(n * p)

reps = 100000

for i = 1:reps
    y = vec(DiffusionSimulation.simulate_data(tdm, taxa))
    global μ .+= y
    global V .+= y * y'
end

V ./= reps
μ ./= reps

μ_tol = 1e-2
v_tol = 5e-2

μ_true = zeros(n * p)
Ψ = RTrees.tree_variance(tree, taxa)
V_true = kron(Σ, Ψ) + kron(Γ, Diagonal(ones(n)))

V_diff = abs.(V - V_true)
μ_diff = abs.(μ - μ_true)

v_max = maximum(V_diff)
μ_max = maximum(μ_diff)

@test v_max < v_tol
@test μ_max < μ_tol
