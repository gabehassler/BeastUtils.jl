module TreeUtils

export rtree,
       standardize_height!

# Internal library for working with trees. functions rely heavily on PhyloNetworks.jl


using PhyloNetworks, Distributions, LinearAlgebra
pn = PhyloNetworks

# for backward compatibility with earlier code
function parse_newick(newick::String)
    return readTopology(newick)
end

@deprecate parse_newick(newick) PhyloNetworks.readTopology(newick)


### Generate HybridTree with random topologies / branch lenghts

# recursively add nodes to HybridNetwork
function random_tips!(net::HybridNetwork,
                        parent::pn.Node,
                        rng::UnitRange{Int}
                        )


    split_ind = rand(first(rng):(last(rng) - 1))

    l_rng = first(rng):split_ind
    r_rng = (split_ind + 1):last(rng)

    num_internal = parent.number

    for new_rng in (l_rng, r_rng)
        if length(new_rng) > 1
            new_parent = pn.Node(num_internal - 1, false)

            pn.parseTreeNode!(new_parent, parent, net)

            random_tips!(net, new_parent, new_rng)

            num_internal -= 2 * length(new_rng) - 1
        else
            leaf_node = pn.Node(new_rng[1], true)
            pn.parseTreeNode!(leaf_node, parent, net)
        end
    end

end

function rtree(n::Int;
                keep_order::Bool = false, # keep (true) or scramble (false) initial label ordering
                ultrametric::Bool = false, # all leaf nodes at same time (true) or not (false)
                edge_dist::ContinuousUnivariateDistribution = Exponential(1.0)
                )
    labels = ["t$i" for i = 1:n]
    return rtree(labels, keep_order = keep_order, ultrametric = ultrametric,
                    edge_dist = edge_dist)
end

# generates a random HybridNetwork
function rtree(labels::Array{<:AbstractString}; # leaf labels
                keep_order::Bool = false, # keep (true) or scramble (false) initial label ordering
                ultrametric::Bool = false, # all leaf nodes at same time (true) or not (false)
                edge_dist::ContinuousUnivariateDistribution = Exponential(1.0)) # distribution for generating random edge lengths
                                                                                # (only applies with ultrametric = false)

    n = length(labels)


    if !ultrametric

        net = HybridNetwork()

        root = pn.Node(-2, false)


        random_tips!(net, root, 1:n)

        pn.pushNode!(net, root)
        net.root = 2 * n - 1

        for edge in net.edge
            edge.length = rand(edge_dist)
        end


    else


        # code largely copied from nj.jl in PhyloNetworks.jl

        available_nodes = [pn.Node(i, true) for i = 1:n]

        node_heights = zeros(n)
        t = 0.0

        net = HybridNetwork(copy(available_nodes), pn.Edge[])
        for i = 1:(n - 1)
            m = n + 1 - i
            i1, i2 = sample(1:m, 2, replace = false)

            n1 = available_nodes[i1]
            n2 = available_nodes[i2]
            t += rand(Exponential(1/m)) #new maximum height

            t1 = t - node_heights[i1]
            t2 = t - node_heights[i2]

            e1 = pn.Edge(net.numEdges + 1, t1)
            e2 = pn.Edge(net.numEdges + 2, t2)

            new_node = pn.Node(net.numNodes + 1, false, false, [e1, e2])

            pn.setNode!(e1, pn.Node[n1, new_node])
            pn.setNode!(e2, pn.Node[n2, new_node])

            pn.setEdge!(n1, e1)
            pn.setEdge!(n2, e2)

            pn.pushEdge!(net, e1)
            pn.pushEdge!(net, e2)

            pn.pushNode!(net, new_node)

            available_nodes[i1] = new_node
            node_heights[i1] = t

            deleteat!(available_nodes, i2) #TODO: make more memory efficient with @view
            deleteat!(node_heights, i2)
        end

        net.root = net.numNodes

    end

    directEdges!(net) # not sure this is necessary
    preorder!(net) # not sure this is necessary

    tip_order = collect(1:n)

    if !keep_order
        tip_order = sortperm(rand(n))
    end

    for i = 1:n
        nm = labels[tip_order[i]]
        net.leaf[i].name = nm
        push!(net.names, nm)
    end

    return net
end


function traverse_node_distances!(dists::Vector{Float64},
                                net::HybridNetwork,
                                node::PhyloNetworks.Node,
                                next_leaf::Int
                                )

    start_ind = next_leaf
    for edge in node.edge
        child = edge.isChild1 ? edge.node[1] : edge.node[2]
        if child !== node
            if child.leaf
                # @show net.leaf[next_leaf].name
                # @show child.name
                @assert net.leaf[next_leaf] === child
                dists[next_leaf] = edge.length
                next_leaf += 1
            else
                rng = traverse_node_distances!(dists, net, child, next_leaf)
                for i in rng
                    dists[i] += edge.length
                end

                next_leaf = last(rng) + 1
            end
        end
    end
    return start_ind:(next_leaf - 1)
end


function leaf_distances(net::HybridNetwork)
    root = net.node[net.root]
    dists = zeros(length(net.leaf))

    traverse_node_distances!(dists, net, root, 1)
    return dists
end

function vcv(net::HybridNetwork, taxa::Vector{String})
    v = PhyloNetworks.vcv(net)
    sim_taxa = string.(names(v))

    perm = indexin(taxa, sim_taxa)
    V = Matrix(v)[perm, perm]
end

function scale_to!(net::HybridNetwork, height::Float64)
    original_height = maximum(getNodeAges(net))
    mult = height / original_height
    for edge in net.edge
        edge.length *= mult
    end
end

function standardize_height!(net::HybridNetwork)
    scale_to!(net, 1.0)
end


end
