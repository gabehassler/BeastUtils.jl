abstract type LoadingsPriorXMLElement <: MyXMLElement end


mutable struct IntegratedFactorsXMLElement <: ModelExtensionXMLElement
    el::XMLOrNothing
    loadings_el::XMLOrNothing
    loadings_prior_els::AbstractArray{XMLElement}
    precision_prior_el::XMLOrNothing
    loadings::S where S <: AbstractArray{Float64, 2}
    precision::S where S <: AbstractArray{Float64, 1}
    precision_scale::Float64
    precision_shape::Float64
    treeModel::TreeModelXMLElement
    mbd::MBDXMLElement
    tree_trait_ind::Int
    standardize_traits::Bool
    msls::Union{MatrixShrinkageLikelihoods, Nothing}
    id::String



end

function IntegratedFactorsXMLElement(treeModel_el::TreeModelXMLElement,
            mbd_el::MBDXMLElement, k::Int; trait_ind::Int = 1)


    p = treeModel_el.trait_dims[trait_ind]

    L = default_loadings(k, p)

    precision = ones(p)
    precision_scale = 1.0
    precision_shape = 1.0
    return IntegratedFactorsXMLElement(nothing, nothing, XMLElement[], nothing,
                                L,
                                precision, precision_scale,
                                precision_shape,
                                treeModel_el,
                                mbd_el, trait_ind, true, nothing,
                                bn.DEFAULT_IF_NAME)
end

function copy(x::IntegratedFactorsXMLElement)
    return IntegratedFactorsXMLElement(
        x.el,
        x.loadings_el,
        x.loadings_prior_els,
        x.precision_prior_el,
        x.loadings,
        x.precision,
        x.precision_scale,
        x.precision_shape,
        x.treeModel,
        x.mbd,
        x.tree_trait_ind,
        x.standardize_traits,
        x.msls,
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


function make_xml(ifxml::IntegratedFactorsXMLElement;
                    reference_precision::Bool = false)
    ifxml.loadings_el = make_loadings(ifxml.loadings)
    if !isnothing(ifxml.msls)
        ifxml.loadings_prior_els = make_xml(ifxml.msls, ifxml.loadings_el)
    else
        ifxml.loadings_prior_els = [make_loadings_normal_prior(ifxml)]
    end


    el = new_element(bn.INTEGRATED_FACTORS)
    attrs = [(bn.ID, ifxml.id),
            (bn.TRAIT_NAME, ifxml.treeModel.node_traits[ifxml.tree_trait_ind])]
    set_attributes(el, attrs)

    if ifxml.standardize_traits
        set_attribute(el, bn.STANDARDIZE, bn.TRUE)
    end

    l_el = new_child(el, bn.LOADINGS)
    add_ref_el(l_el, ifxml.loadings_el)

    p_el = new_child(el, bn.PRECISION)
    if reference_precision
        add_ref_el(p_el, bn.PARAMETER, bn.FACTOR_PRECISION)
    else
        add_parameter(p_el, id=bn.FACTOR_PRECISION, value = ifxml.precision,
                    lower="0")
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
    return ifxml.loadings_el
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
    add_ref_el(d_el, ifxml.loadings_el)
    dist_el = new_child(el, bn.DISTRIBUTION)

    ndm_el = new_child(dist_el, bn.NORMAL_DISTRIBUTION_MODEL)

    m_el = new_child(ndm_el, bn.MEAN)
    add_parameter(m_el, value = [0.0])
    std_el = new_child(ndm_el, bn.STDEV)
    add_parameter(std_el, value =[1.0], lower="0")

    return el
end

function get_priors(xml::IntegratedFactorsXMLElement)
    make_xml(xml)
    if isnothing(xml.msls)
        return [xml.loadings_prior_els; xml.precision_prior_el]
    else
        return [get_priors(xml.msls); xml.precision_prior_el]
    end
end

function get_loggables(xml::IntegratedFactorsXMLElement)
    make_xml(xml)

    loggables = [xml.loadings_el, xml.el[bn.PRECISION][1][bn.PARAMETER][1]]

    if !isnothing(xml.msls)
        make_xml(xml.msls, xml.loadings_el)
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
            msls.shapes[2:end] .= shapes
        end
        if length(scales) > 0
            msls.scales[2:end] .= scales
        end
    end
end
