mutable struct NormalGammaPrecisionOperatorXMLElement <: OperatorXMLElement
    el::XMLOrNothing
    ext_xml::ModelExtensionXMLElement
    td_xml::TraitLikelihoodXMLElement
    weight::Float64

    function NormalGammaPrecisionOperatorXMLElement(
                ext_xml::ModelExtensionXMLElement,
                td_xml::TraitLikelihoodXMLElement)
        return new(nothing, ext_xml, td_xml, 1.0)
    end
end

function make_xml(ngpxml::NormalGammaPrecisionOperatorXMLElement)
    make_xml(ngpxml.ext_xml)
    make_xml(ngpxml.td_xml)

    el = new_element(bn.NGP_OPERATOR)
    set_attribute(el, bn.WEIGHT, ngpxml.weight)
    prior_el = new_child(el, bn.PRIOR)
    add_ref_el(prior_el, get_precision_prior(ngpxml.ext_xml))

    ne_el = new_child(el, bn.NORMAL_EXTENSION)
    set_attribute(ne_el, bn.TREE_TRAIT_NAME, ngpxml.td_xml.attrs[bn.TRAIT_NAME])
    add_ref_el(ne_el, ngpxml.ext_xml.el)
    add_ref_el(ne_el, ngpxml.td_xml.el)

    ngpxml.el = el
    return el
end
