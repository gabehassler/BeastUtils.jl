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
    trait_name::String
    tp_name::String
    standardize_traits::Bool
    msls::Union{MatrixShrinkageLikelihoods, Nothing}


    function IntegratedFactorsXMLElement(treeModel_el::TreeModelXMLElement,
                mbd_el::MBDXMLElement, k::Int)


        p = treeModel_el.trait_dims[1]

        L = default_loadings(k, p)

        trait_name = treeModel_el.node_traits[1]
        tp_name = treeModel_el.param_names[1]
        precision = ones(p)
        precision_scale = 1.0
        precision_shape = 1.0
        tree_trait_ind = 1
        return new(nothing, nothing, XMLElement[], nothing,
                L,
                precision, precision_scale,
                precision_shape,
                treeModel_el,
                mbd_el, tree_trait_ind, trait_name, tp_name, true, nothing)
    end
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


function make_xml(ifxml::IntegratedFactorsXMLElement)
    ifxml.loadings_el = make_loadings(ifxml.loadings)
    if !isnothing(ifxml.msls)
        ifxml.loadings_prior_els = make_xml(ifxml.msls, ifxml.loadings_el)
    else
        ifxml.loadings_prior_els = [make_loadings_normal_prior(ifxml)]
    end


    el = new_element(bn.INTEGRATED_FACTORS)
    attrs = [(bn.ID, bn.DEFAULT_IF_NAME),
            (bn.TRAIT_NAME, ifxml.trait_name)]
    set_attributes(el, attrs)

    if ifxml.standardize_traits
        set_attribute(el, bn.STANDARDIZE, bn.TRUE)
    end

    l_el = new_child(el, bn.LOADINGS)
    add_ref_el(l_el, ifxml.loadings_el)

    p_el = new_child(el, bn.PRECISION)
    add_parameter(p_el, id=bn.FACTOR_PRECISION, value = ifxml.precision,
                    lower="0")
    add_ref_el(el, ifxml.treeModel.el)
    tp_el = new_child(el, bn.TRAIT_PARAMETER)
    add_ref_el(tp_el, bn.PARAMETER, ifxml.tp_name)



    ifxml.el = el

    ifxml.precision_prior_el = make_precision_prior(ifxml)


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
