mutable struct RepeatedMeasuresXMLElement <: ModelExtensionXMLElement
    el::XMLOrNothing
    prior_el::XMLOrNothing
    precision::S where S <: AbstractArray{Float64, 2}
    prior_scale::T where T <: AbstractArray{Float64, 2}
    treeModel_el::TreeModelXMLElement
    mbd_el::MBDXMLElement
    tree_trait_ind::Int
    trait_name::String
    tp_name::String
    standardize_traits::Bool


    function RepeatedMeasuresXMLElement(treeModel_el::TreeModelXMLElement,
                mbd_el::MBDXMLElement)

        p = treeModel_el.trait_dims[1]
        trait_name = treeModel_el.node_traits[1]
        tp_name = treeModel_el.param_names[1]
        precision = Diagonal(ones(p))
        prior_scale = Diagonal(ones(p))
        tree_trait_ind = 1
        return new(nothing, nothing, precision, prior_scale,
                treeModel_el, mbd_el, tree_trait_ind, trait_name, tp_name, true)
    end
end

function make_xml(rm_el::RepeatedMeasuresXMLElement)

    # repeatedMeasures xml element
    el = new_element(bn.REPEATED_MEASURES)
    attrs = [(bn.ID, bn.DEFAULT_RM_NAME), (bn.TRAIT_NAME, rm_el.trait_name)]
    if rm_el.standardize_traits
        set_attribute(el, bn.STANDARDIZE, bn.TRUE)
    end
    set_attributes(el, attrs)
    treeModel_el = make_xml(rm_el.treeModel_el)
    add_ref_el(el, treeModel_el)
    tp_el = new_child(el, bn.TRAIT_PARAMETER)
    add_ref_el(tp_el, bn.PARAMETER, rm_el.tp_name)
    sp_el = new_child(el, bn.SAMPLING_PRECISION)
    mp_el = add_matrix_parameter(sp_el, rm_el.precision, id = bn.DEFAULT_RM_PREC_NAME)

    make_xml(rm_el.mbd_el)

    add_ref_el(el, rm_el.mbd_el.el)
    rm_el.el = el

    # prior xml element

    rm_el.prior_el = make_Wishart_prior(rm_el.prior_scale, mp_el, bn.DEFAULT_RM_PREC_PRIOR_NAME)

    return el
end

function get_loggables(el::RepeatedMeasuresXMLElement)
    make_xml(el)
    return [el.el[bn.SAMPLING_PRECISION][1][bn.MATRIX_PARAMETER][1]]
end

function get_matrix_parameter(el::RepeatedMeasuresXMLElement)
    return get_loggable(el)
end

function get_priors(xml::RepeatedMeasuresXMLElement)
    make_xml(xml)
    return [xml.prior_el]
end

function get_precision_prior(xml::RepeatedMeasuresXMLElement)
    return xml.prior_el
end
