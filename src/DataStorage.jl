module DataStorage

export TraitData,
       csv_to_traitdata,
       parse_traitdata,
       traitdata_to_df

using DataFrames, CSV, BeastUtils.MatrixUtils


struct TraitData
    taxa::Vector{<:AbstractString}
    trait_names::Vector{<:AbstractString}
    data::Matrix{Float64}

    function TraitData(taxa::Vector{<:AbstractString},
                        trait_names::Vector{<:AbstractString},
                        data::Matrix{Float64})
        n, p = size(data)
        if length(taxa) != n
            error("supplied $(length(taxa)) taxa but data has $n rows.")
        end
        if length(trait_names) != p
            error("supplied $(length(trait_names)) taxa but data has $p rows.")
        end
        return new(taxa, trait_names, data)
    end
end

function make_sparse!(X::Matrix{Float64}, sparsity::Float64)

    n, p = size(X)
    p_buffer = zeros(p)

    for i = 1:n

        all_missing = true
        p_buffer .= X[i, :]

        for j = 1:p
            if rand() < sparsity
                X[i, j] = NaN
            else
                all_missing = false
            end
        end

        if all_missing #don't want a taxon with no observations
            ind = rand(1:p)
            X[i, ind] = p_buffer[ind]
        end
    end
end

function store_data(path::String, taxa::Vector{String}, data::Matrix{Float64})
    df = data_to_df(taxa, data)
    CSV.write(path, df)
end

function store_parameters(paths::Array{String, N},
                            params::Array{Array{T, M}, N} where {T, M}) where N

    @assert size(paths) == size(params)

    n = length(paths)


    for i = 1:n
        write(paths[i], join(params[i], ' '))
    end
end

function data_to_df(taxa::Vector{String}, data::Matrix{Float64})
    col_names = ["trait_$i" for i = 1:size(data, 2)]
    return data_to_df(taxa, data, col_names)
end

function data_to_df(taxa::Vector{String}, data::Matrix{Float64}, col_names::Vector{String})
    n, p = size(data)
    @assert length(taxa) == n
    @assert length(col_names) == p
    df = DataFrame()
    df.taxon = taxa

    mis_inds = MatrixUtils.findnans(data)
    mis_data = convert(Matrix{Union{Float64, Missing}}, data)
    mis_data[mis_inds] .= missing

    for i = 1:p
        df[!, Symbol(col_names[i])] = mis_data[:, i]
    end

    return df
end

function df_to_traitdata(df::DataFrame)
    nms = names(df)
    nm1 = lowercase(string(nms[1]))
    if !(nm1 == "taxon" || nm1 == "traits")
        error("The first column name must be 'taxon' or 'traits'.")
    end
    # @assert string(nms[1]) == "taxon"
    n, p = size(df)

    p = p - 1

    data = zeros(n, p)
    for j = 1:p
        for i = 1:n
            val = df[i, j + 1]
            if ismissing(val)
                data[i, j] = NaN
            else
                data[i, j] = val
            end
        end
    end
    td = TraitData(convert(Vector{String}, df[!, 1]), nms[2:end], data)
    df = nothing # there's a memory leak somewhere that this might help with
    return td
end

function df_to_data(df::DataFrame)
    td = df_to_traitdata(df)
    return td.taxa, dt.data
end

function csv_to_data(path::String)
    td = csv_to_traitdata(path)
    return td.taxa, td.data
end

function traitdata_to_df(td::TraitData)
    df = DataFrame(taxon = td.taxa)
    for i = 1:length(td.trait_names)
        df[!, td.trait_names[i]] = td.data[:, i]
    end
    return df
end

function csv_to_traitdata(path::String)
    df = DataFrame(CSV.File(path, missingstrings=["", "NA"]))
    return df_to_traitdata(df)
end

function parse_traitdata(path::String)
    delims = Dict("csv" => ',', "tsv" => '\t', "txt" => '\t')
    ext = path[(findlast(isequal('.'), path) + 1):end]
    delim = delims[ext]
    df = CSV.read(path, DataFrame, delim=delim, missingstring=["", "NA"])
    for i = 2:size(df, 2) #not sure why CSV isn't properly parsing some columns
        if eltype(df[!, i]) <: AbstractString
            df[!, i] = parse.(Float64, df[!, i])
        end
    end
    td = df_to_traitdata(df)
    df = nothing
    return td
end


function data_to_csv(path::String, taxa::Vector{String}, data::Matrix{Float64}, col_names::Vector{String})
    df = data_to_df(taxa, data, col_names)
    CSV.write(path, df)
end

function traitdata_to_csv(path::String, td::TraitData)
    data_to_csv(path, td.taxa, td.data, td.trait_names)
end


end
