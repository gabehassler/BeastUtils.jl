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

mutable struct FactorValidationXMLElement <: MyXMLElement
    el::XMLOrNothing
    ifl::IntegratedFactorsXMLElement
    tdl::TraitLikelihoodXMLElement
    id::String
end

function FactorValidationXMLElement(ifl::IntegratedFactorsXMLElement,
                                    tdl::TraitLikelihoodXMLElement
                                    )
    return FactorValidationXMLElement(nothing, ifl, tdl, "factorValidation")
end

mutable struct CrossValidationXMLElement <: MyXMLElement
    el::XMLOrNothing
    validationProvider_el::MyXMLElement
    log_sum::Bool
    validation_type::String

    CrossValidationXMLElement(xml::MyXMLElement) =
            new(nothing, xml, false, bn.SQUARED_ERROR)
end

function set_validation_type!(cv::CrossValidationXMLElement, t::String)
    cv.validation_type = t
end


function make_xml(cv_el::CrossValidationXMLElement)

    el = new_element(bn.CROSS_VALIDATION)
    set_attributes(el, [(bn.ID, cv_el.validationProvider_el.id),
                        (bn.LOG_SUM, string(cv_el.log_sum)),
                        (bn.TYPE, cv_el.validation_type)])

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

function make_xml(fac_val::FactorValidationXMLElement)
    el = new_element(bn.FACTOR_VALIDATION)
    make_xml(fac_val.ifl)
    make_xml(fac_val.tdl)
    add_ref_el(el, fac_val.ifl.el)
    add_ref_el(el, fac_val.tdl.el)
    trait_name = attribute(fac_val.tdl.el, bn.TRAIT_NAME)
    set_attribute(el, bn.TRAIT_NAME, trait_name)
    fac_val.el = el
end


function get_loggables(cv_el::CrossValidationXMLElement)
    make_xml(cv_el)
    return [cv_el.el]
end
