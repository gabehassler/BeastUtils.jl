module XMLConstructor

export make_MBD_XML


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

include(joinpath(element_dir, "MatrixParameter.jl"))
include(joinpath(element_dir, "DataXMLElement.jl"))
include(joinpath(element_dir, "NewickXMLElement.jl"))
include(joinpath(element_dir, "TreeModelXMLElement.jl"))
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





function get_id(el::XMLElement)
    return attribute(el, bn.ID)
end

function set_id!(el::XMLElement, id::String)
    set_attribute(el, bn.ID, id)
end

function make_Wishart_prior(scale::AbstractArray{Float64, 2},
            mp_el::XMLElement,
            id::String)

    el = new_element(bn.MULTIVARIATE_WISHART_PRIOR)
    set_attributes(el, [(bn.ID, id), (bn.DF, string(size(scale, 1)))])
    scale_el = new_child(el, bn.SCALE_MATRIX)
    add_matrix_parameter(scale_el, scale)
    data_el = new_child(el, bn.DATA)
    add_ref_el(data_el, mp_el)
    return el
end

function add_matrix_parameter(pel::XMLElement, M::AbstractArray{T, 2};
        id::String = "") where T <: Number
    mat_el = new_element(bn.MATRIX_PARAMETER)
    if id != ""
        set_attribute(mat_el, bn.ID, id)
    end
    n, p = size(M)
    for i = 1:n
        param_el = new_element(bn.PARAMETER)
        set_attribute(param_el, bn.VALUE, join(M[i, :], ' '))
        add_child(mat_el, param_el)
    end
    add_child(pel, mat_el)
    return mat_el
end

function add_diagonal_matrix(pel::XMLElement, M::Vector{T};
        id::String = "", lower::String = "") where T <: Number

    mat_el = new_child(pel, bn.DIAGONAL_MATRIX)
    p_el = add_parameter(mat_el, value = M, id=id)
    if lower != ""
        set_attribute(p_el, bn.LOWER, lower)
    end
    return mat_el
end


function add_ref_el(pel::XMLElement, el::XMLElement;
            new_name::String = name(el))
    ref_el = reference_element(el, new_name)
    add_child(pel, ref_el)
    return ref_el
end

function add_ref_el(pel::XMLElement, name::String, id::String)
    ref_el = new_child(pel, name)
    set_attribute(ref_el, bn.IDREF, id)
end

function add_ref_els(pel::XMLElement, els::Array{XMLElement})
    for el in els
        add_ref_el(pel, el)
    end
end

function reference_element(el::XMLElement)
    nm = name(el)
    return reference_element(el, nm)
end

function reference_element(el::XMLElement, nm::String)
    id = attribute(el, bn.ID, required = true)
    ref_el = new_element(nm)
    set_attribute(ref_el, bn.IDREF, id)
    return ref_el
end

function xml_vec(n::Int)
    v = Vector{XMLOrNothing}(undef, n)
    fill!(v, nothing)
    return v
end


struct VarPropXMLElement
    el::XMLOrNothing
    #TODO: more
end


## TODO: move elsewhere

function set_screen_logEvery(bx::BEASTXMLElement, sle::Int)
    mcmc_el = get_mcmc(bx)
    set_screen_logEvery(mcmc_el, sle)
end

function set_file_logEvery(bx::BEASTXMLElement, fle::Int)
    mcmc_el = get_mcmc(bx)
    set_file_logEvery(mcmc_el, fle)
end

function set_filename(bx::BEASTXMLElement, filename::AbstractString)
    mcmc_el = get_mcmc(bx)
    set_filename(mcmc_el, string(filename))
end

function set_data_dates(bx::BEASTXMLElement, dates::AbstractVector{Float64})
    data_el = get_data(bx)
    set_dates(data_el, dates)
end

function get_data(bx::BEASTXMLElement)
    return find_element(bx, DataXMLElement)
end

function get_newick(bx::BEASTXMLElement)
    return find_element(bx, NewickXMLElement)
end

function get_treeModel(bx::BEASTXMLElement)
    return find_element(bx, TreeModelXMLElement)
end

function get_mbd(bx::BEASTXMLElement)
    return find_element(bx, MBDXMLElement)
end

function get_extension(bx::BEASTXMLElement)
    return find_element(bx, ModelExtensionXMLElement)
end

function get_repeatedMeasures(bx::BEASTXMLElement)
    return find_element(bx, RepeatedMeasuresXMLElement)
end

function get_integratedFactorModel(bx::BEASTXMLElement)
    return find_element(bx, IntegratedFactorsXMLElement)
end

function get_traitLikelihood(bx::BEASTXMLElement)
    return find_element(bx, TraitLikelihoodXMLElement)
end

function get_operators(bx::BEASTXMLElement)
    return bx.operators_el.els
end

function get_mcmc(bx::BEASTXMLElement)
    return find_element(bx, MCMCXMLElement)
end

function add_loggable(bx::BEASTXMLElement, el::MyXMLElement;
                        already_made::Bool = false)

    if isnothing(bx.mcmc_el)
        error("MCMC element not created. Can't add loggables")
    end

    if isnothing(bx.loggables)
        bx.loggables = LoggablesXMLElement()
    end

    add_loggable(bx.mcmc_el.loggables, el, already_made = already_made)
    add_loggable(bx.loggables, el, already_made = already_made)
end

function set_full_eval(bx::BEASTXMLElement, n_eval::Int)
    mcmc = get_mcmc(bx)
    mcmc.attrs[bn.FULL_EVALUATION] = string(n_eval)
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
    add_child(bx.el, mbd_el.prior_el)
end

function add_el(bx::BEASTXMLElement, rm_el::RepeatedMeasuresXMLElement)
    make_xml(rm_el)
    add_child(bx.el, rm_el.el)
    add_child(bx.el, rm_el.prior_el)
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
end

function make_xml(bx::BEASTXMLElement)
    xdoc = XMLDocument()
    bx.el = new_element(bn.BEAST)
    set_root(xdoc, bx.el)

    for el in bx.components
        add_el(bx, el)
    end

    return xdoc
end

#




end
