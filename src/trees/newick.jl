const NEWICK_ERROR = "Error parsing newick"
const TAXON_REGEX = r"\s*([^:,\)]+)?\s*(:\s*([0-9\.\-]+))?\s*[,\)]"
const ROOT_REGEX = r"\s*([^:]+)?\s*(:\s*([0-9\.\-]+))?\s*;"

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
        tb = recurse_newick(newick, tree, 2, 1, n_tips + 2, n_tips + 1, 1)
        root_label = parse_root(newick, tb.stop)
        tree.node_labels[1] = root_label
    catch e
        @warn "Parser failed"
        msg = sprint(showerror, e, catch_backtrace())
        println(msg)
    end

    return tree
end

function check_empty(newick::T where T <: AbstractString, ind::Int;
                    end_chars::Array{Char} = [',', ')'])
    searching = true
    valid = true
    while searching
        char = newick[ind]
        if char in end_chars
            searching = false
        elseif !isspace(char)
            valid = false
            searching = false
        end
        ind += 1
    end

    return valid


end

function parse_taxon(newick:: T where T <: AbstractString, offset::Int)

    @show newick[offset:end]

    m = match(TAXON_REGEX, newick, offset)

    @show m

    if isnothing(m[1]) && isnothing(m[3])
        valid = check_empthy(newick, offset)
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

function parse_root(newick::T where T<: AbstractString, offset::Int)

    @show newick[offset:end]

    m = match(ROOT_REGEX, newick, offset)
    @show m
    if isnothing(m[1]) && isnothing(m[3])
        valid = check_empty(newick, offset, end_chars=[';'])
        if !valid
            error("$NEWICK_ERROR: Unable to parse taxon at root.")
        end
    end

    if !isnothing(m[3]) && parse(Float64, m[3]) != 0.0
        @warn "Removing illegal root edge length."
    end

    if !isnothing(m[1])
        return m[1]
    end

    return ""

end

function recurse_newick(newick::T where T<:AbstractString, tree::Tree,
                        start::Int,
                        current_tip::Int,
                        current_internal::Int,
                        parent::Int,
                        edge_ind::Int)
    println("in\n")

    @show tree
    @show current_tip
    @show current_internal
    @show parent
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
        if isspace(char) || char == ',' # ignore whitespace and commas (commas already taken care of in regex)
            ind += 1
            continue
        end

        if char == '('

            println("Entering '('")

            this_node = current_internal
            @show this_node
            this_edge = edge_ind

            tree.edges[edge_ind, 1] = parent
            tree.edges[edge_ind, 2] = this_node

            tb = recurse_newick(newick, tree, ind + 1, current_tip, this_node + 1, this_node, edge_ind + 1) # add others

            current_internal += tb.n_internal

            label, branch_length, ind = parse_taxon(newick, tb.stop)
            @show this_node
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
    @show "out"
    println('\n')

    return TreeBuilder(ind + 1, n_tips, n_internal)


end
