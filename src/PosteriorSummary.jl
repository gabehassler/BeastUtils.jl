module PosteriorSummary

export effective_sample_size,
       hpd_interval,
       hpd_excludes_value

import MCMCDiagnostics.effective_sample_size

function hpd_interval(x::AbstractVector{Float64}; alpha::Float64 = 0.05, sorted::Bool = false)
    n = length(x)
    span = Int(round((1.0 - alpha) * n))
    smallest_rng = Inf
    rng_upper = Inf
    rng_lower = -Inf
    sx::AbstractVector{Float64} = sorted ? x : sort(x)
    lower = 0
    for upper = span:n
        lower = upper - span + 1
        rng = sx[upper] - sx[lower]
        if rng < smallest_rng
            smallest_rng = rng
            rng_lower = sx[lower]
            rng_upper = sx[upper]
        end
    end

    return (lower = rng_lower, upper = rng_upper)
end

function hpd_excludes_value(x::AbstractVector{Float64}, value::Float64; alpha::Float64 = 0.05, sorted::Bool = false)
    hpd = hpd_interval(x, alpha = alpha, sorted = sorted)
    return hpd[1] > value || hpd[2] < value
end

end
