module Logs

using Statistics, DelimitedFiles, BeastUtils.MatrixUtils, JLD, CSV, DataFrames

export LazyLog,
       BeastLog,
       get_log,
       get_log_match,
       import_log,
       list_cols,
       make_meanmatrix,
       make_meanvector,
       make_diagonal_vector,
       make_vector,
       make_correlationMatrix,
       confidence_intervals,
       HDP_intervals,
       log2corr,
       find_cols,
       get_cols,
       combine_logs,
       condense_logs,
       condense_trees,
       make_log,
       compress_log,
       compress_log!,
       read_time

const DELIM_CHAR = '\t'


mutable struct LazyLog
    log_path::String
    df::DataFrame
    burnin::Float64

    function LazyLog(path::String)
        return new(path, DataFrame(), NaN)
    end
end

const BeastLog = Union{String, LazyLog}

function lazy_import!(lz::LazyLog, burnin::Float64)
    if isempty(lz.df) || lz.burnin != burnin
        df = import_log(lz.log_path, burnin = burnin)
        lz.df = df
        lz.burnin = burnin
    end
    return nothing
end

function get_log(lz::LazyLog; burnin::Float64 = 0.0)
    lazy_import!(lz, burnin)
    return names(lz.df), Matrix(lz.df)
end

function import_log(lz::LazyLog; burnin::Float64 = 0.0)
    lazy_import!(lz, burnin)
    return lz.df
end


function compress_log(log_path::String, compressed_path::String)
    cols, data = get_log(log_path, burnin = 0.0)
    save(compressed_path, "cols", cols, "data", data)
end

function compress_log!(log_path::String, compressed_path::String)
    compress_log(log_path, compressed_path)
    rm(log_path)
end

function uncompress_log(compressed_path::String, log_path::String)
    d = load(compressed_path)
    cols = d["cols"]
    data = d["data"]
    make_log(log_path, data, cols, includes_states = true)
end

function uncompress_log!(compressed_path::String, log_path::String)
    uncompress_log(compressed_path, log_path)
    rm(compressed_path)
end

function uncompress_log!(path::String)
    split_path = split(path, '.')
    ext = last(split_path)
    if ext != "jld"
        error("Compressed file must have `.jld` file extension.")
    end

    split_path[end] = "log"
    uncompress_log!(path, join(split_path, '.'))
end




function make_log(path::String, data::AbstractMatrix{Float64},
            col_labels::Vector{String}; includes_states::Bool = false)
    n, p = size(data)
    if length(col_labels) != p
        error("The number of column labels must be the same as the number of columns.")
    end
    open(path, "w") do f
        if !includes_states
            write(f, "state$DELIM_CHAR")
        end
        write(f, join(col_labels, DELIM_CHAR))
        write(f, '\n')
        for i = 1:n
            if !includes_states
                write(f, "$(i - 1)$DELIM_CHAR")
            end
            write(f, join(data[i, :], DELIM_CHAR))
            write(f, '\n')
        end
    end
end

# function below should eventually deprecate get_log()
function import_log(path::String; burnin::Float64 = 0.0, verbose::Bool = true)

    if verbose
        print("importing BEAST log file $path with ")
        printstyled("burnin = $burnin\n", color=:yellow)
    end

    header_row = find_header_row(path)
    all_lines = countlines(path)
    n = all_lines - header_row
    skipto = header_row + Int(round(n * burnin)) + 1

    df = CSV.read(path,
                  DataFrame,
                  header=header_row,
                  skipto = skipto,
                  delim='\t')

    return df
end


function get_log(inpath::String; burnin::Float64 = 0.1)
    if endswith(inpath, ".jld")
        d = load(inpath)
        cols = d["cols"]
        data = d["data"]
        n = size(data, 1)
        n_burn = Int(round(burnin * n))

        return cols, data[(n_burn + 1):end, :]
    else

        col_labels, labels_ind = get_cols_and_ind(inpath)
        p = length(col_labels)
        f = open(inpath)
        m = countlines(f) - labels_ind
        close(f)
        n_burn = Int(round(burnin * m))
        n = m - n_burn
        data = readdlm(inpath, skipstart = labels_ind + n_burn, dims = (n, p))
        return string.(col_labels), data
    end
end

function read_last_line(inpath::String)
    col_labels, labels_ind = get_cols_and_ind(inpath)
    p = length(col_labels)
    f = open(inpath)
    m = countlines(f)

    data = readdlm(inpath, skipstart = m-1, dims = (1, p))
    return string.(col_labels), vec(data)
end

function get_log(inpath::String, cols::Array{Int}; burnin::Float64 = 0.1)
    col_labels, data = get_log(inpath, burnin = burnin)
    return col_labels[cols], data[:, cols]
end

function get_log(inpath::String, header::Union{Regex, String};
                    burnin::Float64 = 0.1)
    Base.depwarn("get_log(path::String, header::$(typeof(header))) is deprecated. " *
                "Use get_log_match(path, header) instead.", :get_log)
    all_cols = get_cols(inpath)
    f, l = find_cols(all_cols, header)
    cols, data = get_log(inpath, collect(f:l), burnin = burnin)
end



function get_log_match(inpath::String, header::Union{Regex, String};
                        burnin::Float64 = 0.1)
    all_cols = get_cols(inpath)
    inds = findall(x -> matches(x, header), all_cols)
    cols, data = get_log(inpath, inds, burnin = burnin)
end

function get_cols(inpath::String)
    if endswith(inpath, ".jld")
        cols = load(inpath, "cols")
        return cols
    else
        cols, ind = get_cols_and_ind(inpath)
        return cols
    end
end

function get_cols_and_ind(inpath::String)
    f = open(inpath)
    col_labels = ""
    ind = 0
    for line in eachline(f)
        ind += 1
        if startswith(line, "state")
            col_labels = line
            break
        end
    end
    close(f)
    cols = split(col_labels)
    return cols, ind
end

function find_header_row(inpath::String)
    counter = 0
    open(inpath, "r") do f
        for line in eachline(inpath)
            counter += 1
            if startswith(line, "state")
                break
            end
        end
    end
    return counter
end


function list_cols(cols::Vector{String})
    p = length(cols)
    for i = 1:p
        println("$i: $(cols[i])")
    end
end

function find_cols(cols::Vector{S}, col::Union{String, Regex}) where S <: AbstractString
    first = 0
    last = 0
    p = length(cols)
    for i = 1:p
        if matches(cols[i], col)
            first = i
            break
        end
    end
    if first == 0
        error("$col not found in column names.")
    end
    for i = p:-1:first
        if matches(cols[i], col)
            last = i
            break
        end
    end
    return first, last
end

function matches(col::S, header::String) where S <: AbstractString
    if startswith(col, header)
        return true
    end
    return false
end

function matches(col::S, header::Regex) where S <: AbstractString
    if occursin(header, col)
        return true
    end
    return false
end

function make_meanmatrix(data::Matrix{T}, cols::Vector{String},
        col::String) where T <: Real

    first, last = find_cols(cols, col)
    return make_meanmatrix(data, first, last)
end


function make_meanmatrix(data::Matrix{T}, first::Int, last::Int) where T <: Real
    t_data = data[:, (first:last)]
    means = mean(t_data, dims = 1)
    q = length(means)
    q = Int(sqrt(q))
    mat = reshape(means, q, q)
    return mat
end

function make_meanmatrix(data::Matrix{T}) where T <: Real
    means = mean(data, dims = 1)
    q = length(means)
    q = Int(sqrt(q))
    mat = reshape(means, q, q)
    return mat
end

function make_meanvector(data::Matrix{T}, first::Int, last::Int) where T <: Real

    means = [mean(data[:, i]) for i in first:last]

    return means
end

function make_diagonal_vector(data::Matrix{T}, first::Int, last::Int) where T <: Real
    q = last - first + 1
    p = Int(sqrt(q))
    n = size(data, 1)
    new_mat = zeros(n, p)
    for i = 1:p
        new_mat[:, i] .= data[:, first + (i - 1) * (p + 1)]
    end
    return new_mat
end

function make_vector(data::Matrix{T}, first::Int, last::Int) where T <: Real
    t_data = data[:, (first:last)]
    return t_data
end

function make_correlationMatrix(X::Matrix{T}) where T <: Real
    n, p = size(X)
    @assert n == p
    C = zeros(n, n)
    for i = 1:n
        for j = 1:n
            C[i, j] = X[i, j] / sqrt(X[i, i] * X[j, j])
        end
    end
    return C
end

function confidence_intervals(d::Matrix{Float64}, p::Int)
    x = d[:, p]
    n = length(x)
    sort!(x)
    bottom = Int(floor(0.025 * n))
    top = Int(ceil(0.975 * n))
    return x[bottom], x[top]
end

function HPD_intervals(d::Matrix{Float64}, p::Int; conf::Float64 = 0.95)
    x = d[:, p]
    return HDP_intervals(x, conf = conf)
end

function HPD_intervals(y::Vector{Float64}; conf::Float64 = 0.95)
    n = length(y)
    x = sort(y)
    s = Int(ceil(conf * n))
    smallest = Inf
    smallest_inds = [0, 0]
    for i = 1:(n - s + 1)
        diff = x[i + s - 1] - x[i]
        if diff < smallest
            smallest = diff
            smallest_inds[1] = i
            smallest_inds[2] = i + s - 1
        end
    end
    return x[smallest_inds]
end

function HPD_intervals(d::Matrix{Float64}, ps::Vector{Int}; conf::Float64 = 0.95)
    hdps = zeros(length(ps), 2)
    for i = 1:length(ps)
        hdpi = HDP_intervals(d, ps[i], conf = conf)
        hdps[i, 1] = hdpi[1]
        hdps[i, 2] = hdpi[2]
    end
    return hdps
end

function log2corr(X::Matrix{Float64}; invert::Bool = false)
    n, q = size(X)
    p = Int(sqrt(q))
    M = zeros(p, p)
    new_X = zeros(n, q)
    for i = 1:n
        for j = 1:q
            M[j] = X[i, j]
        end
        if invert
            M = inv(M)
        end
        corrs = cov2corr(M)
        for j = 1:q
            new_X[i, j] = corrs[j]
        end
    end
    return new_X
end

function get_header(path::String)
    f = open(path)
    lines = Vector{String}(undef, 0)
    for line in eachline(f)
        if startswith(line, "state")
            break
        end
        push!(lines, line)
    end
    close(f)
    return lines
end

function combine_logs(paths::Vector{String}, outpath::String;
        force::Bool = false)
    n = length(paths)
    cols = get_cols(paths[1])
    for i = 2:n
        test_cols = get_cols(paths[i])
        @assert test_cols == cols
    end
    lines = get_header(paths[1])
    push!(lines, join(cols, "\t"))
    state = 0
    all_data = Matrix{Float64}(undef, 0, length(cols))
    for i = 1:n
        cols, data = get_log(paths[i], burnin = 0.0)
        m = size(data, 1)
        for i = 1:m
            data[i, 1] += state
        end
        if i != 1
            data = data[2:end, :]
        end
        all_data = [all_data; data]
        state = data[end, 1]
    end
    m = size(all_data, 1)
    string_data = Vector{String}(undef, m)
    for i = 1:m
        string_data[i] = join(all_data[i, :], "\t")
    end
    lines = [lines; string_data]
    lines = join(lines, "\n")
    write(outpath, lines)
end

function writeline(io::IOStream, s::AbstractString)
    write(io, s)
    write(io, "\n")
end

function condense_logfiles(in_path::String, out_path::String;
                       every::Int = 10,
                       is_valid = is_valid_logfile
    )

    counter = 0

    open(in_path, "r") do input
        open(out_path, "w") do output

            for line in eachline(input)
                if is_valid(line)
                    if rem(counter, every) == 0
                        writeline(output, line)
                    end
                    counter += 1
                else
                    writeline(output, line)
                end
            end
        end
    end
end


function is_valid_logfile(line::AbstractString)::Bool
    first_tab = findfirst('\t', line)
    if isnothing(first_tab)
        return false
    end

    sub_line = @view line[1:(first_tab - 1)]

    try parse(Int, sub_line)
        return true
    catch
        return false
    end
end

function is_valid_treefile(line::AbstractString)::Bool
    return startswith(line, "tree STATE_")
end

function condense_logs(in_path::String, out_path::String; every::Int = 10)
    return condense_logfiles(in_path, out_path, every = every, is_valid = is_valid_logfile)
end

function condense_trees(in_path::String, out_path::String; every::Int = 10)
    return condense_logfiles(in_path, out_path, every = every, is_valid = is_valid_treefile)
end



const TIME_DICT = Dict("seconds" => 1.0,
                       "minutes" => 60.0,
                       "hours" => 3600.0,
                       "days" => 86400.0)


function read_time(path::String; unit::String = "seconds")
    t, u = split(read(path, String))
    t = parse(Float64, t)
    return t * TIME_DICT[u] / TIME_DICT[unit]
end





end
