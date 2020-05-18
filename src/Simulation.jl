module Simulation

using PhyloNetworks, LinearAlgebra, DataFrames

mutable struct TraitSimulationModel
    treeModel::TreeDiffusionModel
    extensionModel::Untion{Nothing, ModelExtension}

    function TraitSimulationModel(treeModel::TreeDiffusionModel,
                                    extensionModel::ModelExtension)
        return new(treeModel, extensionModel)
    end

    function TraitSimulationModel(treeModel::TreeDiffusionModel)
        return new(treeModel, nothing)
    end
end

mutable struct TreeDiffusionModel
    tree::HybridNetwork
    Σ::AbstractArray{Float64, 2}
    μ::AbstractArray{Float64, 1}

    function TreeDiffusionModel(tree::HybridNetwork,
                                Σ::AbstractArray{Float64, 2})
        p, q = size(Σ, 1)
        if p != q
            error("Σ must be square.")
        end
        μ = zeros(p)
        return new(tree, Σ, μ)
    end

    function TreeDiffusionModel(newick::String,
                                Σ::AbstractArray{Float64, 2})
        tree = readTopology(newick)
        return TreeDiffusionModel(tree, Σ)
    end

end

function recuse_tree_diffusion!(df::Matrix{Float64}, )

function simulate_on_tree(model::TreeDiffusionModel)
    p = length(model.μ)
    tree = model.tree
    n = tree.numTaxa
    taxa = [leaf.name for leaf in tree.leaf]
    data = zeros(n, p)
    root = tree.nodes[tree.root]

    recurse_tree_diffusion!(data, root, μ)

    return df
end

function simulate_traits(tsm::TraitSimulationModel)

    data = simulate_on_tree(tsm.treeModel)

    if !isnothing(tsm.extensionModel)
        add_extension!(data, tsm.extensionModel)
    end

    return data

end




end
