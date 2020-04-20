mutable struct TraitValidationXMLElement <: MyXMLElement
    el::XMLOrNothing
    treeModel_el::TreeModelXMLElement
    traitLikelihood_el::TraitLikelihoodXMLElement
    mis_ind::Int
    obs_ind::Int
    standardize::Bool
    id::String

    function TraitValidationXMLElement(tm_el::TreeModelXMLElement,
                    tl_el::TraitLikelihoodXMLElement)

        @assert length(tm_el.node_traits) == 2

        mis_ind = 1
        obs_ind = 2

        return new(nothing, tm_el, tl_el, mis_ind, obs_ind, false, bn.CROSS_VALIDATION)
    end
end

mutable struct CrossValidationXMLElement <: MyXMLElement
    el::XMLOrNothing
    validationProvider_el::TraitValidationXMLElement
    log_sum::Bool

    CrossValidationXMLElement(tv_el::TraitValidationXMLElement) =
            new(nothing, tv_el, false)
end


function make_xml(cv_el::CrossValidationXMLElement)

    el = new_element(bn.CROSS_VALIDATION)
    set_attributes(el, [(bn.ID, cv_el.validationProvider_el.id),
                        (bn.LOG_SUM, string(cv_el.log_sum))])

    make_xml(cv_el.validationProvider_el)
    add_child(el, cv_el.validationProvider_el.el)

    cv_el.el = el
    return el
end


function make_xml(tv_el::TraitValidationXMLElement)

    obs_ind = tv_el.obs_ind
    mis_ind = tv_el.mis_ind

    el = new_element(bn.TRAIT_VALIDATION_PROVIDER)
    attrs = [(bn.ID, bn.TRAIT_VALIDATION),
            (bn.TRAIT_NAME, tv_el.treeModel_el.node_traits[obs_ind]),
            (bn.INFERRED_TRAIT, tv_el.treeModel_el.node_traits[mis_ind]),
            (bn.STANDARDIZE, string(tv_el.standardize))]
    set_attributes(el, attrs)

    make_xml(tv_el.traitLikelihood_el)
    add_ref_el(el, tv_el.traitLikelihood_el.el)

    tp_el = new_child(el, bn.TRAIT_PARAMETER)
    add_ref_el(tp_el, bn.PARAMETER, tv_el.treeModel_el.param_names[obs_ind])

    tv_el.el = el
    return el
end

function get_loggable(cv_el::CrossValidationXMLElement)
    make_xml(cv_el)
    return cv_el.el
end
