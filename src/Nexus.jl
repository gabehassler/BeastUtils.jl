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

end