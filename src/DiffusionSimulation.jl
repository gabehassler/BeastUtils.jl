module DiffusionSimulation

using BeastUtils.RTrees, Distributions, LinearAlgebra

AbsArr = AbstractArray{Float64, 2}

struct TreeDiffusionModel
    tree::RTrees.PhyloTree
    Σ::AbsArr #Diffusion variance
    Γ::AbsArr #Residual variance
    μ0::Vector{Float64}

    function TreeDiffusionModel(newick::String, Σ::AbsArr, Γ::AbsArr)
        @assert isposdef(Σ)
        @assert isposdef(Γ)
        @assert size(Σ) == size(Γ)

        μ0 = zeros(size(Σ, 1))

        return new(RTrees.parse_newick(newick), Σ, Γ, μ0)
    end

    function TreeDiffusionModel(tree::RTrees.PhyloTree, Σ::AbsArr, Γ::AbsArr, μ0::Vector{Float64})
        @assert isposdef(Σ)
        @assert isposdef(Γ)
        @assert size(Σ) == size(Γ)

        return new(tree, Σ, Γ, μ0)
    end


end

function simulate_data(tdm::TreeDiffusionModel, taxa::Vector{String})

    @assert length(tdm.tree.tip_labels) == length(taxa)

    taxa_order = arrange_taxa(tdm.tree.tip_labels, taxa)

    p = size(tdm.Σ, 1)

    data = zeros(tdm.tree.n_tips, p)

    simulate_on_tree!(tdm, data, tdm.tree.root, tdm.μ0, taxa_order) #simulate MBD process

    res_mvn = MvNormal(zeros(p), tdm.Γ) #simulate residual process

    data .+= rand(res_mvn, tdm.tree.n_tips)' #add residual error to tree
end

#Recursively simulates MBD process on tree and fills in 'data' at the tips
function simulate_on_tree!(tdm::TreeDiffusionModel, data::Matrix{Float64},
                node::Int, node_value::Vector{Float64}, taxa_order::Vector{Int})

    tree = tdm.tree
    Σ = tdm.Σ
    n = tree.n_tips

    if node <= n
        error("This should only be called on internal nodes")
    end

    child_inds = RTrees.find_children_inds(tree, node)

    for i in child_inds
        edge_length = tree.edge_lengths[i]
        child = tree.edges[i, 2]

        child_value = copy(node_value)

        if edge_length > 0.0
            mvn = MvNormal(node_value, edge_length * Σ)
            child_value .= rand(mvn)
        elseif edge_length < 0.0
            error("Negative branch length")
        end
        if child <= n
            data[taxa_order[child], :] .= child_value
        else
            simulate_on_tree!(tdm, data, child, child_value, taxa_order)
        end
    end
end

function arrange_taxa(start_order::Vector{String}, target::Vector{String})
    n = length(target)
    @assert length(start_order) == n
    order = zeros(Int, n)
    for i = 1:n
        ind = findfirst(x -> x == target[i], start_order)
        order[ind] = i
    end
    return order
end

end
