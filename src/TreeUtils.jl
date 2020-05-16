module TreeUtils

# Internal library for working with trees. functions rely heavily on PhyloNetworks.jl


using PhyloNetworks, Distributions
pn = PhyloNetworks

function parse_newick(newick::String)
    return readTopology(newick)
end

@deprecate parse_newick(newick) PhyloNetworks.readTopology(newick)


function random_tips!(net::HybridNetwork,
                        parent::pn.Node,
                        rng::UnitRange{Int},
                        leaf_labels::Array{String}
                        )


    split_ind = rand(first(rng):(last(rng) - 1))

    l_rng = first(rng):split_ind
    r_rng = (split_ind + 1):last(rng)

    num_internal = parent.number

    for new_rng in (l_rng, r_rng)
        if length(new_rng) > 1
            new_parent = pn.Node(num_internal - 1, false)

            pn.parseTreeNode!(new_parent, parent, net)

            random_tips!(net, new_parent, new_rng, leaf_labels)

            num_internal -= 2 * length(new_rng) - 1
        else
            leaf_node = pn.Node(new_rng[1], true)
            leaf_node.name = leaf_labels[new_rng[1]]
            push!(net.names, leaf_node.name)
            pn.parseTreeNode!(leaf_node, parent, net)
        end
    end

end

function random_edges!(net::HybridNetwork,
                        edge_dist::ContinuousUnivariateDistribution)

    for edge in net.edge
        edge.length = rand(edge_dist)
    end
end


function rtree(n::Int;
                labels::Array{String} = ["t$i" for i = 1:n],
                keep_order::Bool = false,
                ultrametric::Bool = false,
                edge_dist::ContinuousUnivariateDistribution = Exponential(1.0))

    if length(labels) != n
        throw(ArgumentError("The `labels` argument must be an array of length" *
            " n (only tip labels are currently supported)."))
    end


    tip_order = collect(1:n)

    if !keep_order
        tip_order = sortperm(rand(n))
    end


    net = HybridNetwork()

    if !ultrametric

        root = pn.Node(-2, false)


        random_tips!(net, root, 1:n, labels[tip_order])

        pn.pushNode!(net, root)
        net.root = 2 * n - 1

        random_edges!(net, edge_dist)


    else
        # TODO
    end


    return net
end


# function get_root(edges::Matrix{Int})
#     @assert size(edges, 2) == 2
#     n = size(edges, 1)
#     found = fill(false, n + 1)
#
#     for i = 1:n
#         found[edges[i, 2]] = true
#     end
#
#     not_found = findall(x -> !x, found)
#     @assert length(not_found) == 1
#
#     return not_found[1]
# end
#
# function make_newick(tree::PhyloTree)
#     edge = tree.edges
#     edge_length = tree.edge_lengths
#     Nnode = tree.n_internal
#     tip_labels = tree.tip_labels
#     @rput edge
#     @rput edge_length
#     @rput Nnode
#     @rput tip_labels
#     R"""
#     phy <- list(edge = edge)
#     phy$tip.label <- tip_labels
#     phy$edge.length <- edge_length
#     phy$Nnode <- Nnode
#     class(phy) <- "phylo"
#     newick <- write.tree(phy)
#     """
#     @rget newick
#     return newick
# end
#
# function put_phylo(tree::PhyloTree)
#     edge = tree.edges
#     edge_length = tree.edge_lengths
#     Nnode = tree.n_internal
#     tip_labels = tree.tip_labels
#     @rput edge
#     @rput edge_length
#     @rput Nnode
#     @rput tip_labels
#     R"""
#     phy <- phy <- list(edge = edge)
#     phy$tip.label <- tip_labels
#     phy$edge.length <- edge_length
#     phy$Nnode <- Nnode
#     class(phy) <- "phylo"
#     """
# end
#
# function trim_tree!(tree::PhyloTree, nodes::Vector{Int})
#     put_phylo(tree)
#     @rput nodes
#     R"""
#     drop.tip(phy, nodes)
#     """
# end
#
# function trim_tree(tree::PhyloTree, nodes::Vector{Int})
#     new_tree = deepcopy(tree)
#     trim_tree!(new_tree, node)
#     return new_tree
# end
#
# function trim_tree!(tree::PhyloTree, label::String)
#     node = findfirst(x -> x == label, tree.tip_labels)
#     trim_tree!(tree, node)
# end
#
# function trim_to_n(tree::PhyloTree, n::Int)
#     @assert 2 <= n <= tree.n_tips
#     new_tree = deepcopy(tree)
#
#     for i = tree.n_tips:-1:(n + 1)
#         node = rand(1:i)
#         trim_tree!(new_tree, node)
#     end
#     return new_tree
# end
#
#
#
# function find_parent(tree::PhyloTree, node::Int)
#     edges = tree.edges
#     for i = 1:size(edges, 1)
#         if edges[i, 2] ==  node
#             return edges[i, 1], i
#         end
#     end
#     error("Node $node is either not a node or the root node.")
# end
#
# function find_sib(tree::PhyloTree, parent::Int, node::Int)
#     edges = tree.edges
#     for i = 1:size(edges, 1)
#         if edges[i, 1] == parent && edges[i, 2] != node
#             return edges[i, 2], i
#         end
#     end
#     error("Could not find sibling of node $node with parent $parent.")
# end
#
# function find_children(tree::PhyloTree, node::Int)
#
#     child_inds = find_children_inds(tree, node)
#     return tree.edges[child_inds, 2]
#
# end
#
# function find_children_inds(tree::PhyloTree, node::Int)
#     child_inds = zeros(Int, 0)
#     for i = 1:size(tree.edges, 1)
#         if tree.edges[i, 1] == node
#             push!(child_inds, i)
#         end
#     end
#     return child_inds
# end
#
# #Finds the distance from any node to the root node
# function distance_to_root(tree::PhyloTree, node::Int)
#
#     if node == tree.root
#         return 0.0
#     end
#
#     parent, pn_ind = find_parent(tree, node)
#     distance = tree.edge_lengths[pn_ind]
#
#     while parent != tree.root
#         parent, pn_ind = find_parent(tree, parent)
#         distance += tree.edge_lengths[pn_ind]
#     end
#
#     return distance
# end
#
# #Finds the distance from any taxon to the root node.
# function distance_to_root(tree::PhyloTree, taxon::String)
#     node = findfirst(x -> x == taxon, tree.tip_labels)
#
#     if isnothing(node)
#         error("Taxon \"$taxon\" not in tree.")
#     end
#
#     return distance_to_root(tree, node)
# end
#
#
# #Adds a node above the root node
# function add_root!(tree::PhyloTree, edge_length::Float64)
#     root = tree.root
#
#     for i = 1:length(tree.edges)
#         if tree.edges[i] >= root
#             tree.edges[i] += 1
#         end
#     end
#
#     new_N = tree.n_tips + tree.n_internal
#     new_edges = zeros(Int, new_N, 2)
#     new_lengths = zeros(new_N)
#
#     new_edges[1:(root - 1), :] .= tree.edges[1:(root - 1), :]
#     new_edges[root, :] .= [root, root + 1]
#     new_edges[(root + 1):end, :] .= tree.edges[root:end, :]
#     new_lengths[1:(root - 1)] .= tree.edge_lengths[1:(root - 1)]
#     new_lengths[root] = edge_length
#     new_lengths[(root + 1):end] .= tree.edge_lengths[root:end]
#
#
#     tree.edges = new_edges
#     tree.edge_lengths = new_lengths
#     tree.n_internal = tree.n_internal + 1
# end
#
# function rtree(taxa::Vector{String}; ultrametric::Bool = false)
#     seed = Int(floor(rand() * 1e8))
#     return rtree(taxa, seed, ultrametric = ultrametric)
# end
#
# function rtree(taxa::Vector{String}, seed::Int; ultrametric::Bool = false)
#     n = length(taxa)
#     @rput n
#     @rput taxa
#     @rput seed
#     @rput ultrametric
#     R"""
#     library(ape)
#     set.seed(seed)
#     if (ultrametric) {
#         tree <- rcoal(n, tip.label=taxa)
#     } else {
#         tree <- rtree(n, tip.label=taxa)
#     }
#     """
#     @rget tree
#     return rTree_to_phyloTree(tree)
# end
#
# function tree_variance_sums(tree::PhyloTree; standardize_tree::Bool = false)
#     diag_sum, all_sum, n = recurse_tree_sums(tree, tree.root)
#
#     @assert n == tree.n_tips
#
#     if standardize_tree
#         max_height = maximum([distance_to_root(tree, i) for i = 1:n])
#         diag_sum = diag_sum / max_height
#         all_sum = all_sum / max_height
#     end
#
#     return diag_sum, all_sum
# end
#
# function tree_variance_sums(newick::String; standardize_tree::Bool = false)
#     return tree_variance_sums(parse_newick(newick),
#                                 standardize_tree = standardize_tree)
# end
#
# function recurse_tree_sums(tree::PhyloTree, node::Int)
#
#     n = tree.n_tips
#     if node <= n
#         return 0.0, 0.0, 1
#
#     else
#         d_sum, a_sum = 0.0, 0.0
#         n = 0
#         child_inds = find_children_inds(tree, node)
#         for i in child_inds
#             child = tree.edges[i, 2]
#             edge_length = tree.edge_lengths[i]
#
#             cd_sum, ca_sum, cn = recurse_tree_sums(tree, child)
#
#             n += cn
#             d_sum += cd_sum + cn * edge_length
#             a_sum += ca_sum + cn^2 * edge_length
#         end
#
#         return d_sum, a_sum, n
#     end
# end
#
# function tree_variance(tree::PhyloTree, taxa::Vector{String})
#
#     @assert length(taxa) == tree.n_tips
#
#
#     n = tree.n_tips
#
#     taxa_order = zeros(Int, n)
#
#     for i = 1:n
#         taxon = taxa[i]
#         tree_ind = findfirst(x -> x == taxon, tree.tip_labels)
#         taxa_order[tree_ind] = i
#     end
#
#     Ψ = zeros(n, n)
#
#     recurse_tree_variance(tree, tree.root, Ψ, taxa_order)
#
#     return Ψ
# end
#
# function recurse_tree_variance(tree::PhyloTree, node::Int, Ψ::Matrix{Float64}, taxa_order::Vector{Int})
#     child_inds = find_children_inds(tree, node)
#
#     all_sub_nodes = Vector{Int}(undef, 0)
#
#     for ind in child_inds
#         child_node = tree.edges[ind, 2]
#
#         sub_nodes = Vector{Int}(undef, 0)
#         if child_node <= tree.n_tips
#             sub_nodes = [child_node]
#
#             ψi = taxa_order[child_node]
#             Ψ[ψi, ψi] = tree.edge_lengths[ind]
#         else
#
#             sub_nodes = recurse_tree_variance(tree, child_node, Ψ, taxa_order)
#
#             m = length(sub_nodes)
#
#             for i = 1:m
#                 nodei = sub_nodes[i]
#                 ψi = taxa_order[nodei]
#                 Ψ[ψi, ψi] += tree.edge_lengths[ind]
#
#                 for j = (i + 1):m
#                     nodej = sub_nodes[j]
#                     ψj = taxa_order[nodej]
#                     Ψ[ψi, ψj] += tree.edge_lengths[ind]
#                     Ψ[ψj, ψi] = Ψ[ψi, ψj]
#                 end
#             end
#         end
#
#
#         all_sub_nodes = [all_sub_nodes; sub_nodes]
#
#     end
#
#     return all_sub_nodes
# end
#
#
#
#
# function rescale_tree!(tree::PhyloTree, height::Float64)
#     original_height = maximum([distance_to_root(tree, taxon)
#                                     for taxon in 1:tree.n_tips])
#
#     mult = height / original_height
#     tree.edge_lengths .*= mult
# end
#
# function standardize_tree!(tree::PhyloTree)
#     rescale_tree!(tree, 1.0)
# end


end
