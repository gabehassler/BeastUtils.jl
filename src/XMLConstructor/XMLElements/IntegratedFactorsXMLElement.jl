abstract type LoadingsPriorXMLElement <: MyXMLElement end


mutable struct IntegratedFactorsXMLElement <: ModelExtensionXMLElement
    el::XMLOrNothing
    loadings_prior_els::AbstractArray{XMLElement}
    precision_prior_el::XMLOrNothing
    loadings::MyXMLElement
    precision::S where S <: AbstractArray{Float64, 1}
    precision_scale::Float64
    precision_shape::Float64
    treeModel::TreeModelXMLElement
    tree_trait_ind::Int
    standardize_traits::Bool
    msls::Union{MatrixShrinkageLikelihoods, Nothing}
    rotate_prior::Bool
    rotation_el::XMLOrNothing
    id::String
end

function IntegratedFactorsXMLElement(treeModel_el::TreeModelXMLElement,
                                     k::Int; trait_ind::Int = 1,
                                     orthonormal::Bool = false)


    p = treeModel_el.trait_dims[trait_ind]

    L = default_loadings(k, p)
    if orthonormal
        load_param = ScaledOrthogonalMatrix(L, "L", "U", "scale")
    else
        load_param = MatrixParameter(L, "L", ["L$i" for i = 1:size(L, 1)])
    end

    precision = ones(p)
    precision_scale = 1.0
    precision_shape = 1.0
    return IntegratedFactorsXMLElement(nothing, XMLElement[], nothing,
                                load_param,
                                precision, precision_scale,
                                precision_shape,
                                treeModel_el,
                                trait_ind, true, nothing,
                                false, nothing,
                                bn.DEFAULT_IF_NAME)
end

function make_xmlelement(model::IntegratedFactorModel, tm::TreeModelXMLElement;
                         ind::Int = 1)
    ifx = IntegratedFactorsXMLElement(tm, tip_dimension(model), trait_ind = ind)
    set_loadings(ifx, model.L)
    set_precision(ifx, 1 ./ model.λ)
    return ifx
end

function copy(x::IntegratedFactorsXMLElement)
    return IntegratedFactorsXMLElement(
        x.el,
        x.loadings_prior_els,
        x.precision_prior_el,
        x.loadings,
        x.precision,
        x.precision_scale,
        x.precision_shape,
        x.treeModel,
        x.tree_trait_ind,
        x.standardize_traits,
        x.msls,
        x.rotate_prior,
        x.rotation_el,
        x.id
    )
end

function default_loadings(k::Int, p::Int)
    @assert k <= p
    L = zeros(k, p)
    for i = 1:k
        for j = i:p
            L[i, j] = 1.0
        end
    end
    # L = randn(k, p)
    return L
end

function set_loadings!(ifxml::IntegratedFactorsXMLElement, L::S) where S <: AbstractArray{Float64, 2}
    set_mat!(ifxml.loadings, L)
end

function get_loadings_param(ifxml::IntegratedFactorsXMLElement)
    return ifxml.loadings
end

function make_xml(ifxml::IntegratedFactorsXMLElement;
                    reference_precision::Bool = false)

    make_xml(ifxml.loadings)
    # ifxml.loadings_el = make_loadings(ifxml.loadings)
    actual_prior = nothing
    if !isnothing(ifxml.msls)
        ifxml.loadings_prior_els = make_xml(ifxml.msls)
        actual_prior = ifxml.msls.ms_el

    else
        ifxml.loadings_prior_els = [make_loadings_normal_prior(ifxml)]
        actual_prior = ifxml.loadings_prior_els[1]
    end

    if ifxml.rotate_prior
        ifxml.rotation_el = make_rotation_xml(ifxml, actual_prior)
        push!(ifxml.loadings_prior_els, ifxml.rotation_el)
    end

    el = new_element(bn.INTEGRATED_FACTORS)
    attrs = [(bn.ID, ifxml.id),
            (bn.TRAIT_NAME, ifxml.treeModel.node_traits[ifxml.tree_trait_ind])]
    set_attributes(el, attrs)

    if ifxml.standardize_traits
        set_attribute(el, bn.STANDARDIZE, bn.TRUE)
    end

    l_el = new_child(el, bn.LOADINGS)
    add_ref_el(l_el, ifxml.loadings.el)

    p_el = new_child(el, bn.PRECISION)
    if reference_precision
        add_ref_el(p_el, bn.PARAMETER, bn.FACTOR_PRECISION)
    else
        add_parameter(p_el, id=bn.FACTOR_PRECISION, value = ifxml.precision,
                    lower=0.0)
    end
    add_ref_el(el, ifxml.treeModel.el)
    tp_el = new_child(el, bn.TRAIT_PARAMETER)
    tp_id = ifxml.treeModel.param_names[ifxml.tree_trait_ind]
    add_ref_el(tp_el, bn.PARAMETER, tp_id)



    ifxml.el = el

    if !reference_precision
        ifxml.precision_prior_el = make_precision_prior(ifxml)
    end


    return el
end


function make_precision_prior(ifxml::IntegratedFactorsXMLElement)
    el = new_element(bn.GAMMA_PRIOR)
    set_attributes(el, [bn.ID => "$(bn.FACTOR_PRECISION).prior",
                        bn.SCALE => string(ifxml.precision_scale),
                        bn.SHAPE => string(ifxml.precision_shape)])
    add_ref_el(el, ifxml.el[bn.PRECISION][1][bn.PARAMETER][1])
    return el
end


function get_hmc_parameter(ifxml::IntegratedFactorsXMLElement)
    make_xml(ifxml)
    return ifxml.loadings.el
end

function make_loadings(L::AbstractArray{Float64, 2})
    el = new_element(bn.MATRIX_PARAMETER)
    set_attribute(el, bn.ID, bn.DEFAULT_LOADINGS_ID)
    n, p = size(L)
    for i = 1:n
        add_parameter(el, id = "$(bn.DEFAULT_LOADINGS_ID)$i",
                        value = L[i, :])
    end
    return el
end

function get_normal_prior(ifxml::IntegratedFactorsXMLElement)
    if ifxml.rotate_prior
        return ifxml.rotation_el
    end
    if isnothing(ifxml.msls)
        return ifxml.loadings_prior_els[1]
    else
        return ifxml.msls.ms_el
    end
end

function make_loadings_normal_prior(ifxml::IntegratedFactorsXMLElement)
    el = new_element(bn.DISTRIBUTION_LIKELIHOOD)
    set_attribute(el, bn.ID, "$(bn.DEFAULT_LOADINGS_ID).prior")

    d_el = new_child(el, bn.DATA)
    make_xml(ifxml.loadings)
    add_ref_el(d_el, ifxml.loadings.el)
    dist_el = new_child(el, bn.DISTRIBUTION)

    ndm_el = new_child(dist_el, bn.NORMAL_DISTRIBUTION_MODEL)

    m_el = new_child(ndm_el, bn.MEAN)
    add_parameter(m_el, value = [0.0])
    std_el = new_child(ndm_el, bn.STDEV)
    add_parameter(std_el, value =[1.0], lower=0.0)

    return el
end

function make_rotation_xml(ifxml::IntegratedFactorsXMLElement,
                           prior::XMLElement)
    el = new_element(bn.NORMAL_ORTHOGONAL_SUBSPACE)
    set_attribute(el, bn.ID, ifxml.id * ".orthogonalPrior")
    add_ref_el(el, ifxml.loadings.el)
    add_ref_el(el, prior)
    return el
end

function get_priors(xml::IntegratedFactorsXMLElement)
    make_xml(xml)

    if isnothing(xml.msls)
        priors = [xml.loadings_prior_els; xml.precision_prior_el]
    else
        priors = [get_priors(xml.msls); xml.precision_prior_el]
    end

    if xml.rotate_prior
        if isnothing(xml.msls)
            replacement_ind = 1
        else
            replacement_ind = findfirst(x -> x === xml.msls.ms_el, priors)
        end
        priors[replacement_ind] = xml.rotation_el
    end

    return priors
end

function get_loggables(xml::IntegratedFactorsXMLElement)
    make_xml(xml)

    loggables = [xml.loadings.el, xml.el[bn.PRECISION][1][bn.PARAMETER][1]]

    if !isnothing(xml.msls)
        make_xml(xml.msls)
        loggables = [loggables; get_loggables(xml.msls)]
    end

    return loggables
end

function get_precision_prior(xml::IntegratedFactorsXMLElement)
    return xml.precision_prior_el
end

function set_shrinkage_mults!(ifxml::IntegratedFactorsXMLElement;
                        shapes::Vector{Float64} = Float64[],
                        scales::Vector{Float64} = Float64[])
    msls = ifxml.msls
    if isnothing(msls)
        error("No shrinkage prior on the integrated factors model.")
    else
        if length(shapes) > 0
            fill_shrinkage_array!(msls.shapes, shapes)
        end
        if length(scales) > 0
            fill_shrinkage_array!(msls.scales, scales)
        end
    end
end

function fill_shrinkage_array!(to_fill::Vector{Float64},
                               fill_from::Vector{Float64})
    n = length(to_fill)
    m = length(fill_from)
    if n == m
        to_fill .= fill_from
    elseif (n - 1) == m
        to_fill[2:end] .= fill_from
    else
        error("Not implemented for these array dimensions.")
    end
end

function set_loadings(pfa::IntegratedFactorsXMLElement,
                      L::AbstractArray{Float64, 2})
    mat_param = pfa.loadings
    set_mat!(mat_param, L)
end


function set_precision(ifa::IntegratedFactorsXMLElement,
                       mat::AbstractArray{Float64, 2})
    if !isdiag(mat)
        error("Factor model must have diagonal precision")
    end
    set_precision(ifa, diag(mat))
end

function set_precision(ifa::IntegratedFactorsXMLElement,
                                precs::AbstractArray{Float64, 1})
    p = length(precs)
    if length(ifa.precision) != p
        error("Cannot set precision. Incompatible dimensions.")
    end
    for i = 1:p
        if precs[i] <= 0
            error("Precisions must be positive.")
        end
    end
    ifa.precision = precs
end