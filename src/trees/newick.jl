const NEWICK_ERROR = "Error parsing newick"
const TAXON_REGEX = r"\s*([^:,\)]+)?\s*(:\s*([0-9\.\-]+))?\s*[,\)]"

struct TreeBuilder
    # start::Int
    stop::Int
    # current_tip::Int
    # current_internal::Int
    n_tips::Int
    n_internal::Int
end

function initial_check_valid(newick::String)
    s_newick = strip(newick)
    err_message = ""

    n_open = 0
    n_closed = 0
    n_comma = 0

    if s_newick[1] != '('
        error("$NEWICK_ERROR: First non-whitespace character must be '('.")
    elseif s_newick[end] != ';'
        error("$NEWICK_ERROR: Last non-whitespace character must be ';'.")
    end

    n_open = count_chars(s_newick, '(')
    n_close = count_chars(s_newick, ')')

    if n_open != n_close
        error("$NEWICK_ERROR: the number of open '(' and close ')' must be the same.")
    end

    n_comma = count_chars(s_newick, ',')

    n_tips = n_comma + 1
    n_internal = n_open

    n_tips, n_internal
end

function count_chars(s::S where S <:AbstractString, c::Char)
    return count(x -> x == c, s)
end


function parse_newick(newick::String)
    n_tips, n_internal = initial_check_valid(newick)

    tree = Tree(n_tips, n_internal)

    # tb = TreeBuilder(1, 0, 1, n_tips + 1, n_tips, n_internal)
    try
        recurse_newick(newick, tree, 2, 1, n_tips + 2, n_tips + 1, 1)
    catch e
        @warn "Parser failed"
        msg = sprint(showerror, e, catch_backtrace())
        println(msg)
    end

    return tree
end

function parse_taxon(newick:: T where T <: AbstractString, offset::Int)

    @show newick[offset:end]

    m = match(TAXON_REGEX, newick, offset)

    @show m

    if isnothing(m[1]) && isnothing(m[3])
        searching = true
        valid = true
        while searching
            char = newick[ind]
            if char == ',' || char == ')'
                searching = false
            elseif !isspace(char)
                valid = false
                searching = false
            end
        end

        if !valid
            error("$NEWICK_ERROR: Unable to parse taxon $(m.match)")
        end
    end

    label = ""
    branch_length = NaN


    if !isnothing(m[1])
        label = m[1]
    end

    if !isnothing(m[3])
        try
            branch_length = parse(Float64, m[3])
        catch
            error("$NEWICK_ERROR: Unable to parse branch length $(m[3]).")
        end
    end

    return label, branch_length, m.offset + length(m.match) - 1
end

function recurse_newick(newick::T where T<:AbstractString, tree::Tree,
                        start::Int,
                        current_tip::Int,
                        current_internal::Int,
                        parent::Int,
                        edge_ind::Int)


    @show current_tip
    @show current_internal
    @show edge_ind

    started = false
    o_count = 0
    c_count = 0

    searching = true

    stop = 0
    n_tips = 0
    n_internal = 0
    ind = start

    group_taxa = 0
    while searching
        char = newick[ind]
        if isspace(char) || char == ',' # ignore whitespace
            ind += 1
            continue
        end

        if char == '('

            println("Entering '('")

            this_node = current_internal
            this_edge = edge_ind

            tree.edges[edge_ind, 1] = parent
            tree.edges[edge_ind, 2] = this_node
            current_internal += 1

            tb = recurse_newick(newick, tree, ind + 1, current_tip, current_internal, this_node, edge_ind + 1) # add others

            label, branch_length, ind = parse_taxon(newick, tb.stop + 1)
            tree.node_labels[this_node - tree.n_tips] = label
            tree.edge_lengths[this_edge] = branch_length

            edge_ind = edge_ind + tb.n_tips + tb.n_internal + 1

            n_tips += tb.n_tips
            n_internal += tb.n_internal + 1 # need to add the current internal node
            current_internal += n_internal
            current_tip += n_tips



        elseif char == ')'
            println("Entering ')'")
            # TODO: parse tree length (and maybe name)
            searching = false
        else
            println("Entering taxon")
            label, branch_length, ind = parse_taxon(newick, ind)

            tree.edges[edge_ind, 1] = parent
            tree.edges[edge_ind, 2] = current_tip
            tree.edge_lengths[edge_ind] = branch_length
            tree.tip_labels[current_tip] = label
            current_tip += 1
            edge_ind += 1
            n_tips += 1

            @show newick[ind]

        end


    end

    return TreeBuilder(stop, n_tips, n_internal)


end
