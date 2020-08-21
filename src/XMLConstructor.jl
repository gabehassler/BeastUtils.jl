module XMLConstructor

export save_xml,
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
       DataModel


using LightXML, LinearAlgebra, DataFrames, PhyloNetworks
import BeastUtils.BeastNames, BeastUtils.TreeUtils
bn = BeastNames

XMLOrNothing = Union{XMLElement, Nothing}

abstract type MyXMLElement end
abstract type SimpleXMLElement <: MyXMLElement end

mutable struct BEASTXMLElement
    el::XMLOrNothing
    components::Vector{MyXMLElement}

    function BEASTXMLElement()
        return new(nothing, MyXMLElement[])
    end
end

abstract type ModelExtensionXMLElement <: MyXMLElement end
abstract type OperatorXMLElement <: MyXMLElement end

const element_dir = joinpath(@__DIR__, "XMLConstructor", "XMLElements")
const dir = joinpath(@__DIR__, "XMLConstructor")

include(joinpath(dir, "constructor_helpers.jl"))

include(joinpath(element_dir, "NothingXMLElement.jl"))
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

import LightXML: find_element
function find_element(bx::BEASTXMLElement, type::Type)
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

    return bx.components[found[1]]
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
        insert!(bx.components, ind, loggables)
    elseif length(old_loggables) == 1
        join_loggables(old_loggables[1], loggables)
    else
        error("Cannot join with more than one LoggablesXMLElement")
    end
end






function make_xml(n::Nothing)
    # do nothing
end

function add_el(bx::BEASTXMLElement, no::Nothing)
    #do nothing
end

function add_el(bx::BEASTXMLElement, el::SimpleXMLElement)
    make_xml(el)
    add_child(bx.el, el.el)
end

function add_el(bx::BEASTXMLElement, data_el::DataXMLElement)
    make_xml(data_el)
    add_child(bx.el, data_el.el)
end

function add_el(bx::BEASTXMLElement, newick_el::NewickXMLElement)
    make_xml(newick_el)
    add_child(bx.el, newick_el.el)
end

function add_el(bx::BEASTXMLElement, treeModel::TreeModelXMLElement)
    make_xml(treeModel)
    add_child(bx.el, treeModel.el)
end

function add_el(bx::BEASTXMLElement, mbd_el::MBDXMLElement)
    make_xml(mbd_el)
    add_child(bx.el, mbd_el.el)
    add_child(bx.el, mbd_el.precision_prior.el)
end

function add_el(bx::BEASTXMLElement, rm_el::RepeatedMeasuresXMLElement)
    make_xml(rm_el)
    add_child(bx.el, rm_el.el)
    add_child(bx.el, rm_el.precision_prior.el)
end

function add_el(bx::BEASTXMLElement, if_el::IntegratedFactorsXMLElement)
    make_xml(if_el)
    add_child(bx.el, if_el.loadings_el)
    for prior_el = if_el.loadings_prior_els
        add_child(bx.el, prior_el)
    end
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

function add_el(bx::BEASTXMLElement, tl_el::TraitLikelihoodXMLElement)
    make_xml(tl_el)
    add_child(bx.el, tl_el.el)
end

function add_el(bx::BEASTXMLElement, os_el::OperatorsXMLElement)
    make_xml(os_el)
    add_child(bx.el, os_el.el)
end

function add_el(bx::BEASTXMLElement, mc_el::MCMCXMLElement)
    make_xml(mc_el)
    add_child(bx.el, mc_el.el)
end

function add_el(bx::BEASTXMLElement, t_el::TimerXMLElement)
    make_xml(t_el)
    add_child(bx.el, t_el.el)
end

function add_el(bx::BEASTXMLElement, xml::TraitLoggerXMLElement)
    make_xml(xml)
    add_child(bx.el, xml.el)
end

function add_el(bx::BEASTXMLElement, el::FactorLogPredictiveDensity)
    make_xml(el)
    add_child(bx.el, el.fac_el)
    add_child(bx.el, el.like_el)
    add_child(bx.el, el.cl_el)
end

function add_el(bx::BEASTXMLElement, el::CrossValidationXMLElement)
    make_xml(el)
    add_child(bx.el, el.el)
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

function add_el(bx::BEASTXMLElement, mat_inv::MatrixInverseXMLElement)
    make_xml(mat_inv)
    add_child(bx.el, mat_inv.el)
end

function add_el(bx::BEASTXMLElement, corr::CorrelationMatrixXMLElement)
    make_xml(corr)
    add_child(bx.el, corr.el)
end





end
