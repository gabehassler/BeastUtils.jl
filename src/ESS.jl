module ESS

using Statistics

export ess


function ess(X::Matrix{Float64}, cols::Vector{Int})
    p = length(cols)
    esses = zeros(p)
    for i = 1:p
        esses[i] = ess(X, cols[i])
    end
    return esses
end

function ess(X::Matrix{Float64})
    p = size(X, 2)
    return ess(X, collect(1:p))
end

function ess(X::Matrix{Float64}, col::Int)
    n = size(X, 1)
    ρ = 0.0
    for i = 1:(n - 3)
        p = lag_correlation(X, col, i)
        # @show i
        # @show p
        if p > 0.0
            ρ += p
        else
            break
        end
        # ρ += p
    end
    return n / (1 + 2 * ρ)
end

# function lag_correlation(X::Matrix{Float64}, col::Int, lag::Int)
#     n = size(X, 1)
#     x1 = X[1:(n - lag), col]
#     x2 = X[(lag + 1):n, col]
#     @assert length(x1) == length(x2)
#     ρ = cor(x1, x2)
#     return ρ
# end

function lag_correlation(X::Matrix{Float64}, col::Int, lag::Int)
    si = 0.0
    sj = 0.0
    ssi = 0.0
    ssj = 0.0
    sij = 0.0
    n = size(X, 1) - lag
    for i = 1:n
        j = i + lag
        x = X[i, col]
        y = X[j, col]

        si += x
        ssi += x * x
        sj += y
        ssj += y * y
        sij += x * y
    end

    ei = si / n
    ej = sj / n
    ei2 = ssi / n
    ej2 = ssj / n
    eij = sij / n

    sdi = sqrt(ei2 - ei * ei)
    sdj = sqrt(ej2 - ej * ej)
    covij = eij - ei * ej
    return covij / (sdi * sdj)
end





end
