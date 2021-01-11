module XMLConstructor

export BEASTXMLElement,
       MCMCOptions,
       save_xml,
       make_residual_xml,
       make_pfa_xml,
       make_joint_xml,
       set_mbd_precision,
       set_residual_precision,
       set_loadings,
       DiffusionModel,
       ResidualVarianceModel,
       IntegratedFactorModel,
       JointProcessModel,
       DataModel,
       tip_dimension,
       data_dimension,
       set_options!,
       set_chain_length!,
       get_operators,
       set_file_logEvery


using LightXML, LinearAlgebra, DataFrames, PhyloNetworks, UnPack
import BeastUtils.BeastNames, BeastUtils.TreeUtils
using BeastUtils.MatrixUtils
bn = BeastNames

XMLOrNothing = Union{XMLElement, Nothing}

abstract type MyXMLElement end
abstract type SimpleXMLElement <: MyXMLElement end
abstract type AbstractDataXMLElement <: MyXMLElement end

mutable struct BEASTXMLElement
    el::XMLOrNothing
    components::Vector{MyXMLElement}
end

mutable struct MCMCOptions
    chain_length::Int
    file_log_every::Int
    screen_log_every::Int
    likelihood_check_count::Int
end

function MCMCOptions(;
                     chain_length::Int = 10_000,
                     file_log_every::Int = maximum(1, div(chain_length, 1_000)),
                     screen_log_every::Int = file_log_every,
                     likelihood_check_count::Int = -1)

    if likelihood_check_count == -1
        check_count = div(chain_length, 10)
        check_count = maximum(check_count, 100)
        check_count = minimum(check_count, 1000)
    else
        check_count = likelihood_check_count
    end

    return MCMCOptions(chain_length, file_log_every, screen_log_every,
                       likelihood_check_count)
end

import Base.show

function Base.show(io::IO, bx::BEASTXMLElement)
    disp = "$(typeof(bx)):\n"
    disp *= "components:"
    for comp in bx.components
        disp *= "\n\t$(typeof(comp))"
    end
    println(io, disp)
end

function BEASTXMLElement()
    return BEASTXMLElement(nothing, MyXMLElement[])
end

abstract type ModelExtensionXMLElement <: MyXMLElement end
abstract type OperatorXMLElement <: MyXMLElement end

const element_dir = joinpath(@__DIR__, "XMLConstructor", "XMLElements")
const dir = joinpath(@__DIR__, "XMLConstructor")

include(joinpath(dir, "constructor_helpers.jl"))
include(joinpath(dir, "integration.jl"))

include(joinpath(element_dir, "NothingXMLElement.jl"))
include(joinpath(element_dir, "CompoundXMLElement.jl"))
include(joinpath(element_dir, "MatrixParameter.jl"))
include(joinpath(element_dir, "DataXMLElement.jl"))
include(joinpath(element_dir, "NewickXMLElement.jl"))
include(joinpath(element_dir, "TreeModelXMLElement.jl"))
include(joinpath(element_dir, "WishartPriorXMLElement.jl"))
include(joinpath(element_dir, "MBDXMLElement.jl"))
include(joinpath(element_dir, "RepeatedMeasuresXMLElement.jl"))
include(joinpath(element_dir, "MatrixShrinkage.jl"))
include(joinpath(element_dir, "IntegratedFactorsXMLElement.jl"))
include(joinpath(element_dir, "TraitLikelihoodXMLElement.jl"))
include(joinpath(element_dir, "LatentFactorModelXMLElement.jl"))
include(joinpath(element_dir, "LogPredictiveDensity.jl"))
include(joinpath(element_dir, "NormalGammaPrecisionOperatorXMLElement.jl"))
include(joinpath(element_dir, "PrecisionGibbsOperatorXMLElement.jl"))
include(joinpath(element_dir, "OperatorsXMLElement.jl"))
include(joinpath(element_dir, "LoggablesXMLElement.jl"))
include(joinpath(element_dir, "MCMCXMLElement.jl"))
include(joinpath(element_dir, "TimerXMLElement.jl"))
include(joinpath(element_dir, "TraitValidationXMLElement.jl"))
include(joinpath(element_dir, "HMCOperatorXMLElement.jl"))
include(joinpath(element_dir, "LoadingsGibbsOperatorXMLElement.jl"))
include(joinpath(element_dir, "OldLoadingsGibbsOperatorXMLElement.jl"))
include(joinpath(element_dir, "LatentFactorPrecisionOperatorXMLElement.jl"))
include(joinpath(element_dir, "FactorTreeGibbsOperatorXMLElement.jl"))
include(joinpath(element_dir, "TraitLoggerXMLElement.jl"))
include(joinpath(element_dir, "LKJPrecisionXMLElement.jl"))
include(joinpath(element_dir, "NormalMatrixNormLikelihood.jl"))
include(joinpath(element_dir, "ScaledOrthogonalMatrix.jl"))
include(joinpath(element_dir, "MultivariateGammaLikelihood.jl"))
include(joinpath(element_dir, "NormalDistributionLikelihood.jl"))
include(joinpath(element_dir, "MultiplicativeScalePrior.jl"))


include(joinpath(dir, "utils.jl"))
include(joinpath(dir, "constructors.jl"))
include(joinpath(dir, "getters.jl"))
include(joinpath(dir, "setters.jl"))




function make_xml(bx::BEASTXMLElement)
    xdoc = XMLDocument()
    bx.el = new_element(bn.BEAST)
    set_root(xdoc, bx.el)

    for el in bx.components
    add_el(bx, el)
    end

    return xdoc
end


function save_xml(path::String, bx::BEASTXMLElement;
                  change_filename::Bool = true)
    nm = basename(path)
    s = split(nm, '.')
    @assert length(s) == 2
    @assert s[2] == "xml"

    if change_filename
        filename = s[1]
        set_filename(bx, filename)
    end

    xdoc = make_xml(bx)
    save_file(xdoc, path)
    free(xdoc)
    return nothing
end




import LightXML: add_child
function add_child(xel::XMLElement, no::Nothing)
    # do nothing
end

function add_child(bx::BEASTXMLElement, el::MyXMLElement)
    push!(bx.components, el)
end

function add_child(bx::BEASTXMLElement, el::MyXMLElement, ind::Int)
    insert!(bx.components, ind, el)
end

import LightXML: find_element

function find_element(bx::BEASTXMLElement, type::Type)
    return bx.components[find_element_ind(bx, type)]
end

function find_element_ind(bx::BEASTXMLElement, type::Type)
    if !(type <: MyXMLElement)
        error("Can only find objects of type MyXMLElement.")
    end

    found = findall(x -> typeof(x) <: type, bx.components)

    if length(found) == 0
        error("No elements of type $type.")
    elseif length(found) > 1
        error("Multiple elements of type $type. " *
              "Please use function 'find_elements'.")
    end

    return found[1]
end


function find_elements(bx::BEASTXMLElement, type::Type)
    found = findall(x -> typeof(x) <: type, bx.components)
    return bx.components[found]
end

function add_parameter_id(pel::XMLElement, id::String)
    el = new_element(bn.PARAMETER)
    set_attribute(el, bn.ID, id)
    add_child(pel, el)
    return el
end




function add_loggable(bx::BEASTXMLElement, el::MyXMLElement;
                        already_made::Bool = false)

    loggables = LoggablesXMLElement([el], [already_made])
    add_loggables(bx, loggables)
end


function add_loggables(bx::BEASTXMLElement, loggables::LoggablesXMLElement)
    old_loggables = find_elements(bx, LoggablesXMLElement)
    if length(old_loggables) == 0
        mcmc = get_mcmc(bx)
        ind = findfirst(x -> x === mcmc, bx.components)
        mcmc_loggables = mcmc.loggables
        push!(mcmc_loggables, loggables)
        insert!(bx.components, ind, mcmc_loggables)
    elseif length(old_loggables) == 1
        push!(old_loggables[1], loggables)
    else
        error("Cannot join with more than one LoggablesXMLElement")
    end
    return nothing
end






function make_xml(n::Nothing)
    # do nothing
end

function add_el(bx::BEASTXMLElement, no::Nothing)
    #do nothing
end

function add_el(bx::BEASTXMLElement, el::TextXMLElement)
    add_child(bx.el, el.el)
end

function add_el(bx::BEASTXMLElement, el::MyXMLElement)
    make_xml(el)
    add_child(bx.el, el.el)
end


function add_el(bx::BEASTXMLElement, mbd_el::MBDXMLElement)
    make_xml(mbd_el)
    add_child(bx.el, mbd_el.el)
    add_el(bx, mbd_el.precision_prior)
end


function add_el(bx::BEASTXMLElement, n_el::NothingXMLElement)
    # do nothing
end

function add_el(bx::BEASTXMLElement, cp_el::CompoundXMLElement)
    for el in cp_el.els
        add_el(bx, el)
    end
end

function add_el(bx::BEASTXMLElement, rm_el::RepeatedMeasuresXMLElement)
    make_xml(rm_el)
    add_child(bx.el, rm_el.el)
    add_child(bx.el, rm_el.precision_prior.el)
end

function add_el(bx::BEASTXMLElement, if_el::IntegratedFactorsXMLElement)
    make_xml(if_el)
    add_el(bx, if_el.loadings)
    add_el(bx, if_el.loadings_prior)
    add_child(bx.el, if_el.el)
    add_child(bx.el, if_el.precision_prior_el)
end

function add_el(bx::BEASTXMLElement, lfm_el::LatentFactorModelXMLElement)
    make_xml(lfm_el)
    add_child(bx.el, lfm_el.loadings_el)
    add_child(bx.el, lfm_el.loadings_prior_el)
    add_child(bx.el, lfm_el.el)
    add_child(bx.el, lfm_el.precision_prior_el)
end

function add_el(bx::BEASTXMLElement, el::FactorLogPredictiveDensity)
    make_xml(el)
    add_child(bx.el, el.fac_el)
    add_child(bx.el, el.like_el)
    add_child(bx.el, el.cl_el)
end

function add_el(bx::BEASTXMLElement, els::Array{MyXMLElement})
    for el in els
        add_el(bx, el)
    end
end

function add_el(bx::BEASTXMLElement, loggables::LoggablesXMLElement)
    for i = 1:length(loggables.els)
        if !loggables.already_made[i]
            l = loggables.els[i]
            add_el(bx, l)
        end
    end
end

function add_el(bx::BEASTXMLElement, el::ScaledOrthogonalMatrix)
    make_xml(el)
    add_child(bx.el, el.scale.el)
    add_child(bx.el, el.U.el)
    add_child(bx.el, el.el)
end

function add_el(bx::BEASTXMLElement, el::MatrixShrinkageLikelihoods)
    make_xml(el)
    for mult in el.mult_els
        add_child(bx.el, mult)
    end
    for gp in el.gp_els
        add_child(bx.el, gp)
    end
    for ls in el.local_scale_els
        add_child(bx.el, ls)
    end
    for gs in el.global_scale_els
        add_child(bx.el, gs)
    end
    for bb in el.bb_els
        add_child(bx.el, bb)
    end
    add_child(bx.el, el.ms_el)
    for gp in el.global_prior_els
        add_child(bx.el, gp)
    end
end

function add_el(bx::BEASTXMLElement, el::MultiplicativeScalePrior)
    add_el(bx, el.mults)
    add_el(bx, el.precs)
    add_el(bx, el.mult_prior)
    add_el(bx, el.scale_prior)
end


end
