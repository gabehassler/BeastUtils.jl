mutable struct OperatorsXMLElement
    el::XMLOrNothing
    els::Vector{OperatorXMLElement}

    OperatorsXMLElement(el::OperatorXMLElement) = new(nothing, [el])
    OperatorsXMLElement(els::Vector{OperatorXMLElement}) = new(nothing, els)

end

function make_xml(os_el::OperatorsXMLElement)
    for el in os_el.els
        make_xml(el)
    end

    el = new_element(bn.OPERATORS)
    set_attribute(el, bn.ID, bn.OPERATORS)
    for o_el in os_el.els
        add_child(el, o_el.el)
    end
    os_el.el = el
    return el
end

mutable struct FireParameterOperatorXMLElement <: OperatorXMLElement
    el::XMLOrNothing
    param_provider::MyXMLElement
    value::Array{T} where T <: Real
    weight::Float64

    FireParameterOperatorXMLElement(my_xml::MyXMLElement) =
            new(nothing, my_xml, zeros(0), 1.0)
end

function make_xml(fpo_el::FireParameterOperatorXMLElement)

    make_xml(fpo_el.param_provider)

    el = new_element(bn.FIRE_PARAMETER_CHANGED)
    set_attribute(el, bn.WEIGHT, fpo_el.weight)

    if length(fpo_el.value) > 0
        value_string = join(fpo_el.value, ' ')
        set_attribute(el, bn.VALUE, value_string)
    end

    param_el = get_matrix_parameter(fpo_el.param_provider)
    add_ref_el(el, param_el)

    fpo_el.el = el

    return el
end
