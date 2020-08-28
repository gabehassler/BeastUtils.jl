module SlowLikelihoods

export TreeData,
       PFAParameters,
       loglikelihood,
       conditional_meanvar,
       tree_var,
       var

using PhyloNetworks, LinearAlgebra, Distributions, UnPack, Statistics
using BeastUtils.MatrixUtils, BeastUtils.Simulation, BeastUtils.XMLConstructor

const DEFAULT_PSS = 0.001

struct TreeData
    taxa::Vector{String}
    data::Matrix{Float64}

    function TreeData(taxa::Vector{String}, data::Matrix{Float64})
        if length(taxa) != size(data, 1)
            error("Incompatible dimensions.")
        end
        return new(taxa, data)
    end
end

mutable struct ModelParameters
    Σ::AbstractMatrix{Float64} # diffusion variance
    Γ::AbstractMatrix{Float64} # residulal variance
    L::AbstractMatrix{Float64} # loadings
    pss::Float64
    μ::AbstractVector{Float64} # prior mean
end

function PFAParameters(L::Matrix{Float64}, γ::Vector{Float64})
    k, p = size(L)
    return ModelParameters(Diagonal(ones(k)),
                           Diagonal(γ),
                           L,
                           Inf,
                           zeros(k)
                          )
end

function tree_var(tree::HybridNetwork, taxa::Vector{String})
    Ψ_df = vcv(tree)
    tree_taxa = names(Ψ_df)
    perm = indexin(taxa, tree_taxa)
    @assert taxa == tree_taxa[perm]

    Ψ = Matrix(Ψ_df)[perm, perm]
    return Ψ
end

import Statistics: var

function var(params::ModelParameters, tree::HybridNetwork, taxa::Vector{String};
             rescale_tree::Bool = false)
    @unpack L, Σ, Γ = params
    Ψ = tree_var(tree, taxa)
    if rescale_tree
        scale_factor = maximum(diag(Ψ))
        Ψ ./= scale_factor
    end
    Ψ .+= 1.0 / params.pss

    n = size(Ψ, 1)

    V = kron(L' * params.Σ * L, Ψ) + kron(params.Γ, Diagonal(ones(n)))

    return V
end

function conditional_meanvar(params::ModelParameters, tree::HybridNetwork,
                             data::TreeData, inds::Vector{Int})
    V = var(params, tree, data.taxa)
    μ = mean(params, data)
    vdata = vec(data.data)
    obs_inds = findall(!isnan, vdata)
    cond_inds = setdiff(obs_inds, inds)
    Mi = @view V[inds, inds]
    Mic = @view V[inds, cond_inds]
    Mcc = V[cond_inds, cond_inds] # can't invert view
    Pcc = inv(Mcc)
    V_cond = Mi - Mic * Pcc * Mic'

    μc = @view μ[cond_inds]
    μi = @view μ[inds]

    yc = @view vdata[cond_inds]

    μ_cond = μi + Mic * (Pcc * (yc - μc))

    return μ_cond, V_cond
end



function mean(params::ModelParameters, data::TreeData)
    n = length(data.taxa)
    k, p = size(params.L)
    μ_fac = zeros(n, k)
    for i = 1:n
        μ_fac[i, :] .= params.μ
    end
    μ = vec(μ_fac * params.L)

    return μ
end

function loglikelihood(params::ModelParameters, tree::HybridNetwork,
                       data::TreeData;
                       rescale_tree::Bool = false)

    V = var(params, tree, data.taxa, rescale_tree = rescale_tree)
    make_symmetric!(V)
    μ = mean(params, data)

    v_data = vec(data.data)
    obs_inds = findall(!isnan, v_data)
    V_obs = @view V[obs_inds, obs_inds]
    data_obs = @view v_data[obs_inds]
    μ_obs = @view μ[obs_inds]

    dist = MvNormal(Vector(μ_obs), Symmetric(V_obs)) # TODO: make memory efficient (but not actually)
    return logpdf(dist, data_obs)
end

function loglikelihood(tsm::TraitSimulationModel, data::Matrix{Float64};
                       rescale_tree::Bool = true)
    params = sim_to_params(tsm)
    tree = tsm.treeModel.tree
    return loglikelihood(params, tree, TreeData(tsm.taxa, data),
                         rescale_tree = rescale_tree)
end

function loglikelihood(jpm::JointProcessModel, dm::DataModel, newick::String;
                       rescale_tree::Bool = true)
    params = xml_to_params(jpm)
    tree = readTopology(newick)
    td = TreeData(dm.taxa, dm.data)
    return loglikelihood(params, tree, td, rescale_tree = rescale_tree)
end

function sim_to_params(tsm::TraitSimulationModel; pss::Float64 = DEFAULT_PSS)
    @unpack tree, Σ, μ = tsm.treeModel
    ext = tsm.extensionModel
    if isnothing(ext)
        p = size(Σ, 1)
        return ModelParameters(Σ, Diagonal(zeros(p)), Diagonal(ones(p)), pss, μ)
    elseif typeof(ext) == LatentFactorModel
        @unpack L, Λ = ext
        return ModelParameters(Σ, Λ, L, pss, μ)
    elseif typeof(ext) == ResidualVarianceModel
        p = size(Σ, 1)
        @unpack Γ = ext
        return ModelParameters(Σ, Γ, Diagonal(ones(p)), pss, μ)
    else
        error("Unkown extension type: $(typeof(ext))")
    end
end

function xml_to_params(jpm::JointProcessModel)
    @unpack μ, Σ, pss = jpm.diffusion_model
    exts = jpm.extensions
    k = tip_dimension(jpm)
    p = data_dimension(jpm)
    L = zeros(k, p)
    Γ = zeros(p, p)

    p_offset = 0
    k_offset = 0
    for ext in exts
        p_ext = data_dimension(ext)
        k_ext = tip_dimension(ext)
        p_rng = (1 + p_offset):(p_ext + p_offset)
        k_rng = (1 + k_offset):(k_ext + k_offset)
        L_sub, Γ_sub = xmlext_to_params(ext)
        L[k_rng, p_rng] .= L_sub
        Γ[p_rng, p_rng].= Γ_sub

        p_offset += p_ext
        k_offset += k_ext
    end

    return ModelParameters(Σ, Γ, L, pss, μ)
end



function xmlext_to_params(ifm::IntegratedFactorModel)
    return ifm.L, Diagonal(ifm.λ)
end

function xmlext_to_params(rvm::XMLConstructor.ResidualVarianceModel)
    p = data_dimension(rvm)
    return Diagonal(ones(p)), rvm.Γ
end




end