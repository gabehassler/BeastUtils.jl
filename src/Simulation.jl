module Simulation

using PhyloNetworks, LinearAlgebra, DataFrames

abstract type ModelExtension end

mutable struct TreeDiffusionModel
    tree::HybridNetwork
    Σ::AbstractArray{Float64, 2}
    μ::AbstractArray{Float64, 1}

    function TreeDiffusionModel(tree::HybridNetwork,
                                Σ::AbstractArray{Float64, 2})
        p, q = size(Σ)
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

mutable struct TraitSimulationModel
    taxa::AbstractArray{AbstractString, 1}
    treeModel::TreeDiffusionModel
    extensionModel::Union{Nothing, ModelExtension}

    function TraitSimulationModel(treeModel::TreeDiffusionModel,
                                    extensionModel::ModelExtension)
        return new(treeModel, extensionModel)
    end

    function TraitSimulationModel(treeModel::TreeDiffusionModel)
        return new(treeModel, nothing)
    end
end



function get_dimension(tdm::TreeDiffusionModel)
    return length(tdm.μ)
end

function simulate_on_tree(model::TreeDiffusionModel,
                        taxa::AbstractArray{T, 1}) where T <: AbstractString
    params = PhyloNetworks.ParamsMultiBM(model.μ, model.Σ)
    trait_sim = simulate(model.tree, params)
    sim_taxa = trait_sim.M.tipNames
    perm = indexin(taxa, sim_taxa)
    @assert sim_taxa[perm] == taxa #TODO remove

    n = length(taxa)
    p = get_dimension(model)

    sim_Y = trait_sim[:Tips]
    Y = sim_Y'[perm, :]

    return Y
end


function simulate(tsm::TraitSimulationModel)

    data = simulate_on_tree(tsm.treeModel)

    if !isnothing(tsm.extensionModel)
        add_extension!(data, tsm.extensionModel)
    end

    return data

end




end
