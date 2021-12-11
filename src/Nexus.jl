module Nexus

export parse_taxon_block

const NTAXA_DIM_START = "Dimensions ntax="

function parse_taxon_block(path::String)

    taxa = String[]

    is_ntaxa_line = false
    current_taxon = 1

    for line in eachline(path)
        stripped_line = strip(line)
        if isempty(line)
            continue
        elseif stripped_line == "Begin taxa;"
            is_ntaxa_line = true
        elseif is_ntaxa_line
            @assert startswith(stripped_line, NTAXA_DIM_START) && endswith(line, ';')
            str_ntaxa = stripped_line[(length(NTAXA_DIM_START) + 1):(end - 1)]
            ntaxa = parse(Int, str_ntaxa)
            taxa = fill("", ntaxa)
            is_ntaxa_line = false
        elseif !isempty(taxa)
            if stripped_line == "Taxlabels"
                continue
            elseif length(taxa) >= current_taxon
                taxa[current_taxon] = stripped_line
                current_taxon += 1
            else # all taxa done
                break;
            end
        end
    end

    return taxa
end

function parse_taxa_dict(path::String)
    taxa = Dict{String, String}()

    past_translate_line = false

    for line in eachline(path)
        s_line = strip(line)
        if past_translate_line && s_line == ";"
            break
        elseif past_translate_line && !isempty(s_line)
            vals = split(s_line)
            taxon_number = vals[1]

            taxon_name = join(vals[2:end], ' ')
            if endswith(taxon_name, ',')
                taxon_name = taxon_name[1:(end - 1)]
            end
            taxa[taxon_name] = taxon_number
        elseif s_line == "Translate"
            past_translate_line = true
        end
    end
    return taxa
end


function findsecond(c::Char, s::AbstractString)
    first_found = false
    for i = 1:length(s)
        if s[i] == c
            if first_found
                return i
            else
                first_found = true
            end
        end
    end

    return nothing
end

import Base.strip
function strip(::Nothing)
    # do nothing
end

struct TreeIterator
    lineIterator::Base.EachLine{IOStream}
end

TreeIterator(path::AbstractString) = TreeIterator(eachline(path))

import Base.iterate

function iterate(iterator::TreeIterator)
    lineIterator = iterator.lineIterator
    next_line = iterate(lineIterator)

    while !isnothing(next_line) && !startswith(next_line[1], "tree STATE_")
        next_line = iterate(lineIterator)
    end

    if isnothing(next_line)
        return nothing
    end

    next_line = next_line[1]

    @assert endswith(next_line, ';')
    state_start = findfirst('_', next_line) + 1
    state_end = findsecond(' ', next_line) - 1
    state = parse(Int, @view next_line[state_start:state_end])
    tree_start = findfirst('(', next_line)
    tree = next_line[tree_start:end] # ends with ; by assertion
    return ((tree, state), nothing)
end

function iterate(iterator::TreeIterator, ::Int)
    iterate(iterator)
end

function iterate(iterator::TreeIterator, ::Nothing)
    iterate(iterator)
end



end