module TreeUtils

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

# generates a random HybridNetwork
function rtree(n::Int; # number of taxa
                labels::Array{String} = ["t$i" for i = 1:n], # leaf labels (automatically [t1, ..., tn])
                keep_order::Bool = false, # keep (true) or scramble (false) initial label ordering
                ultrametric::Bool = false, # all leaf nodes at same time (true) or not (false)
                edge_dist::ContinuousUnivariateDistribution = Exponential(1.0)) # distribution for generating random edge lengths
                                                                                # (only applies with ultrametric = false)

    if length(labels) != n
        throw(ArgumentError("The `labels` argument must be an array of length" *
            " n (only tip labels are currently supported)."))
    end


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
            @show m
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


mutable struct ParamsMultiBM <: pn.ParamsProcess
    mu::AbstractArray{Float64, 1}
    sigma::AbstractArray{Float64, 2}
    randomRoot::Bool
    varRoot::AbstractArray{Float64, 2}
    shift::Union{ShiftNet, Missing}
    L::LowerTriangular{Float64}
end

ParamsMultiBM(mu::AbstractArray{Float64, 1},
                sigma::AbstractArray{Float64, 2}) =
        ParamsMultiBM(mu, sigma, false, Diagonal([NaN]), missing, cholesky(sigma).L)



function process_dim(params::ParamsMultiBM)
    return length(params.mu)
end

function anyShift(params::ParamsMultiBM)
    # TODO
    @warn "Shifts not currently implemented. Returning false."
    return false
end



function initSimulateMBD(nodes::Vector{pn.Node}, params::Tuple{ParamsMultiBM})
    n = length(nodes)
    p = process_dim(params[1])
    vals = zeros(p, n) # random values
    means = zeros(p, n) # means
    return (means, vals)
end

function updateRootSimulateMBD!(X::Tuple{Matrix{Float64}, Matrix{Float64}},
                                i::Int,
                                params::Tuple{ParamsMultiBM})
    params = params[1]
    means = X[1]
    vals = X[2]
    if (params.randomRoot)
        means[:, i] .= params.mu # expectation
        vals[:, i] .= params.mu + params.L * randn(length(params.mu)) # random value #TODO: make memory efficient
    else
        means[:, i] .= params.mu # expectation
        vals[:, i] .= params.mu # random value
    end
end

# Going down to a tree node
function updateTreeSimulateBM!(X::Tuple{Matrix{Float64}, Matrix{Float64}},
                               i::Int,
                               parentIndex::Int,
                               edge::pn.Edge,
                               params::Tuple{ParamsBM})
    params = params[1]
    means = X[1]
    vals = X[2]
    means[:, i] .= means[:, parentIndex] + params.shift.shift[i] # expectation
    vals[:, i] .= vals[:, parentIndex] + params.shift.shift[i] + sqrt(edge.length) * parasm.L * randn(length(params.mu)) # random value #TODO: make memory efficient
end

# Going down to an hybrid node
function updateHybridSimulateBM!(X::Tuple{Matrix{Float64}, Matrix{Float64}},
                                 i::Int,
                                 parentIndex1::Int,
                                 parentIndex2::Int,
                                 edge1::pn.Edge,
                                 edge2::pn.Edge,
                                 params::Tuple{ParamsBM})

    params = params[1]
    means = X[1]
    vals = X[2]
    p = process_dim(params)
    means[:, i] .= edge1.gamma * means[:, parentIndex1] + edge2.gamma * means[:, parentIndex2] # expectation
    vals[:, i] .=  edge1.gamma * (vals[:, parentIndex1] + sqrt(edge1.length) * params.L * randn(p)) +
                    edge2.gamma * (vals[:, parentIndex2] + sqrt(edge2.length) * params.L * randn(p)) # random value
    # TODO: shifts?
    # TODO: memory efficincy
end



end
