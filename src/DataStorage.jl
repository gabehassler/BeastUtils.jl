module DataStorage

using DataFrames, CSV, BeastUtils.MatrixUtils

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

function df_to_data(df::DataFrame)
    nms = names(df)
    @assert string(nms[1]) == "taxon"
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
    return convert(Vector{String}, df[!, 1]), data
end

function csv_to_data(path::String)
    df = CSV.read(path)
    return df_to_data(df)
end


end
