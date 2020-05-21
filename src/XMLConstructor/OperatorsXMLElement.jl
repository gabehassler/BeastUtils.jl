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
        if typeof(o_el) <: CompoundOperatorXMLElement
            for xml_el in o_el.els
                add_child(el, xml_el)
            end
        else
            add_child(el, o_el.el)
        end
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

abstract type CompoundOperatorXMLElement <: OperatorXMLElement end

mutable struct ShrinkageScaleOperators <: CompoundOperatorXMLElement
    els::Vector{XMLOrNothing}
    msl::MatrixShrinkageLikelihoods
    ifxml::IntegratedFactorsXMLElement
    scale_factors::Vector{Float64}
    scale_weights::Vector{Float64}
    bb_weights::Vector{Float64}
    fix_first::Bool

    function ShrinkageScaleOperators(msl::MatrixShrinkageLikelihoods,
                                    ifxml::IntegratedFactorsXMLElement)
        k = get_fac_dim(msl)
        return new(xml_vec(2 * k),
                    msl,
                    ifxml,
                    fill(0.5, k),
                    ones(k),
                    ones(k),
                    false)
    end
end

function make_scale_operator(param::XMLElement,
                            scale_factor::Float64,
                            weight::Float64)
    el = new_element(bn.SCALE_OPERATOR)
    set_attribute(el, bn.SCALE_FACTOR, scale_factor)
    set_attribute(el, bn.WEIGHT, weight)
    add_ref_el(el, param)
    return el
end

function make_bb_operator(bb::XMLElement, weight::Float64)
    el = new_element(bn.BAYESIAN_BRIDGE_GIBBS_OP)
    set_attribute(el, bn.WEIGHT, weight)
    add_ref_el(el, bb)
    return el
end

function make_xml(xml::ShrinkageScaleOperators)
    msl = xml.msl
    make_xml(xml.ifxml)
    make_xml(msl, xml.ifxml.loadings_el)
    k = get_fac_dim(msl)

    if xml.fix_first
        xml.els[1] = nothing
    else
        xml.els[1] = make_scale_operator(msl.gp_els[1], xml.scale_factors[1],
                                        xml.scale_weights[1])
    end

    for i = 2:k
        xml.els[i] = make_scale_operator(msl.mult_els[i - 1],
                                        xml.scale_factors[i],
                                        xml.scale_weights[i])
    end

    for i = 1:k
        xml.els[i + k] = make_bb_operator(msl.bb_els[i], xml.bb_weights[i])
    end

end
