module MatrixUtils

using LinearAlgebra, Statistics, DataFrames

export cov2corr,
    make_symmetric!,
    maxind,
    eye,
    eye_diag,
    matrix_inds,
    make_upper_triangular!,
    pivot_cov,
    pivot_mean,
    findnans,
    replace_nans!,
    permutation_matrix,
    standardize_data!,
    standardize!,
    log!,
    issquare,
    randpd,
    missing_var,
    missing_cov,
    missing_mean,
    missing_cor

function cov2corr(Σ::Matrix{Float64})
    n, p = size(Σ)
    @assert n == p
    P = zeros(n, n)
    for i = 1:n
        for j = i:n
            P[i, j] = Σ[i, j] / sqrt(Σ[i, i] * Σ[j, j])
            P[j, i] = P[i, j]
        end
    end
    return P
end

function cov2corr(Σ::Diagonal{Float64})
    p = size(Σ, 1)
    return Diagonal(ones(p))
end


function make_symmetric!(X::Matrix{T}; use_upper::Bool = true) where T <: Any
    n, m = size(X)
    @assert n == m
    for i = 1:(n - 1)
        for j = (i + 1):n
            if use_upper
                X[j, i] = X[i, j]
            else
                X[i, j] = X[j, i]
            end
        end
    end
end

function maxind(X::AbstractArray{<:Real, 1})
    n = length(X)
    m = 1
    Xm = X[m]
    for i = 2:n
        if X[i] > Xm
            m = i
            Xm = X[m]
        end
    end
    return m
end

function maxind(X::AbstractArray{<:Real, 2})
    ind = maxind(vec(X))
    m, n = size(X)
    col, r = divrem(ind - 1, m)
    col += 1
    row = mod1(r, m)
    return (row, col)
end

function eye(n::Int)
    return diagm(0 => ones(n))
end

function eye_diag(n::Int)
    return Diagonal(ones(n))
end

function matrix_inds(start::Int, stop::Int; diagonals::Bool = true)

    n = stop - start + 1
    p = Int(round(sqrt(n)))
    @assert p^2 == n
    q = div(p * (p + 1), 2)
    if !diagonals
        q = q - p
    end
    inds = zeros(Int, q)
    counter = 0
    for i = 1:p
        r = i:p
        if !diagonals
            r = (i + 1):p
        end
        for j in r
            counter += 1
            v_dim = (i - 1) * p + j
            inds[counter] = v_dim + start - 1
        end
    end
    return inds
end

function matrix_inds(n::Int; diagonals::Bool = true)
    return matrix_inds(1, n^2, diagonals = diagonals)
end

function make_upper_triangular!(X::Matrix{Float64};
        positive_diagonals::Bool = false)
    n, p = size(X)
    for i = 2:n
        for j = 1:(i - 1)
            X[i, j] = 0.0
        end
    end
    if positive_diagonals
        for i = 1:min(n, p)
            X[i, i] = abs(X[i, i])
        end
    end
end

function pivot_cov(M::Hermitian{Float64}, dims::Vector{Int};
        obs_dims::Vector{Int} = setdiff(1:size(M, 1), dims))

    Mmm = M[dims, dims]
    Mmo = M[dims, obs_dims]
    Moo = M[obs_dims, obs_dims]
    return Mmm - Mmo * inv(Moo) * Mmo'
end

function pivot_mean(Y::AbstractArray{Float64, 1}, M::Hermitian{Float64},
        μ::Vector{Float64}, dims::Vector{Int};
        obs_dims::Vector{Int} = setdiff(1:size(M, 1), dims))

    μm = μ[dims]
    Mmo = M[dims, obs_dims]
    Moo = M[obs_dims, obs_dims]
    deviation = Y[obs_dims] - μ[obs_dims]
    return μm + Mmo * (inv(Moo) * deviation)
end

function findnans(x::AbstractArray{Float64}; reverse::Bool = false)
    if reverse
        return findall(y -> !isnan(y), x)
    else
        return findall(y -> isnan(y), x)
    end
end

function permutation_matrix(x::Vector{Int})
    n = length(x)
    Q = zeros(n, n)
    for i = 1:n
        Q[i, x[i]] = 1.0
    end
    return Q
end

function replace_nans!(x::AbstractArray{Float64}; replacement::Float64 = 0.0)
    n = length(x)
    for i = 1:n
        if isnan(x[i])
            x[i] = replacement
        end
    end
end

function standardize_data!(data::Matrix{Float64})
    n, p = size(data)
    ns = zeros(Int, p)
    sums = zeros(p)
    sum_squares = zeros(p)

    for i = 1:n
        for j = 1:p
            x = data[i, j]
            if !isnan(x)
                ns[j] += 1
                sums[j] += x
                sum_squares[j] += x^2
            end
        end
    end

    means = sums ./ ns
    sds = sqrt.(sum_squares ./ ns - means.^2)

    for j = 1:p
        data[:, j] .-= means[j]
        data[:, j] ./= sds[j]
    end
end

function standardize!(x::Vector{Float64})
    μ = 0.0
    σ2 = 0.0
    n_obs = 0
    for i = 1:length(x)
        if !isnan(x[i])
            μ += x[i]
            σ2 += x[i] * x[i]
            n_obs += 1
        end
    end
    if n_obs < 2
        error("Unable to standardize. The number of non-missing elements is $n_obs.")
    end
    μ /= n_obs
    σ2 /= n_obs
    sd = sqrt(σ2 - μ * μ)
    for i = 1:length(x)
        x[i] = (x[i] - μ) / sd
    end
    return x
end

function log!(x::T) where T <: AbstractArray{Float64}
    for i = 1:length(x)
        x[i] = log(x[i])
    end
    return x
end

function log_standardize!(x::T) where T <: AbstractArray{Float64}
    return standardize!(log!(x))
end

function issquare(x::AbstractArray{T, 2}) where T <: Any
    n, m = size(x)
    return n == m
end

function missing_mean(x::AbstractMatrix{Union{Missing, Float64}})
    missing_mean(convert_to_float(x))
end

function missing_mean(x::AbstractMatrix{Float64})
    n, p = size(x)
    μ = zeros(p)
    ns = zeros(Int, p)
    for j = 1:p
        for i = 1:n
            if !isnan(x[i, j])
                μ[j] += x[i, j]
                ns[j] += 1
            end
        end
    end

    return μ ./ ns
end


function missing_var(x::AbstractArray{Union{Float64, Missing}})
    return missing_var(convert_to_float(x))
end

function missing_var(x::AbstractArray{Float64})
    μ = 0.0
    σ2 = 0.0

    n = length(x)
    for i = 1:n
        if !isnan(x[i])
            μ += x[i]
            σ2 += x[i] * x[i]
        end
    end

    μ /= n
    σ2 /= n
    σ2 -= μ * μ

    return σ2
end

function missing_cov(x::AbstractArray{Union{Float64, Missing}})
    return missing_cov(convert_to_float(x))
end

function missing_cov(X::Matrix{Float64})
    n, p = size(X)
    counts = zeros(Int, p, p)
    Σ = zeros(p, p)
    μ = zeros(p, p)
    V = zeros(p, p)
    # Σ_row = zeros(p, p)
    # Σ_col = zeros(p, p)
    # μ = zeros(p)
    # μ_row = zeros(p, p)
    # μ_col = zeros(p, p)
    for i = 1:n
        for j = 1:p
            xj = X[i, j]
            if !isnan(xj)
                Σ[j, j] += xj * xj
                μ[j, j] += xj
                counts[j, j] += 1

                for k = (j + 1):p
                    xk = X[i, k]
                    if !isnan(xk)
                        Σ[j, k] += xj * xk
                        V[j, k] += xj * xj
                        V[k, j] += xk * xk
                        μ[j, k] += xj
                        μ[k, j] += xk
                        counts[j, k] += 1
                    end
                end
            end
        end
    end

    make_symmetric!(counts)
    μ ./= counts
    Σ ./= counts
    V ./= counts

    for i = 1:p
        Σ[i, i] -= μ[i, i] * μ[i, i]
        for j = (i + 1):p
            Σ[i, j] -= μ[i, j] * μ[j, i]
            Σ[j, i] = Σ[i, j]
            V[i, j] -= μ[i, j] * μ[i, j]
            V[j, i] -= μ[j, i] * μ[j, i]
        end
    end


    return Σ, V
end

function missing_cor(X::Matrix{Float64})
    Σ, V = missing_cov(X)
    p = size(Σ, 1)
    for i = 1:p
        Σ[i, i] = 1.0
        for j = (i + 1):p
            Σ[i, j] /= sqrt(V[i, j] * V[j, i])
            Σ[j, i] = Σ[i, j]
        end
    end
    return Σ
end

function missing_cor(x::AbstractArray{<:Union{Float64, Missing}},
                     y::AbstractArray{<:Union{Float64, Missing}})
    x_inds = findall(a -> !ismissing(a) && !isnan(a), x)
    y_inds = findall(a -> !ismissing(a) && !isnan(a), y)
    both_inds = intersect(x_inds, y_inds)
    return cor((@view x[both_inds]), (@view y[both_inds]))
end




function convert_to_float(X::Matrix{Float64})
    return X
end

function convert_to_float(X::Matrix{Union{Missing, Float64}})
    s = size(X)
    Y = zeros(s)
    for i = 1:length(X)
        if ismissing(X[i])
            Y[i] = NaN
        else
            Y[i] = X[i]
        end
    end
    return Y
end


function missing_cor(df::DataFrame)
    M = Matrix(df)
    M = convert_to_float(M)
    return missing_cor(M)
end

function randpd(n::Int; rescale::Bool = true)
    X = randn(n, n)
    M = X' * X
    if rescale
        diags = diag(M)
        v = var(diags, corrected=false)
        sd = sqrt(v)
        M ./= sd
    end
    return Symmetric(M)
end

end
