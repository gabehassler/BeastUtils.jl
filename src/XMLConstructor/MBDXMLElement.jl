mutable struct MBDXMLElement <: MyXMLElement
    el::XMLOrNothing
    prior_el::XMLOrNothing
    precision::AbstractArray{Float64, 2}
    prior_scale::AbstractArray{Float64, 2}
    is_random::Bool
    diagonal_prec::Bool


    MBDXMLElement(p::Int) = new(nothing, nothing, Diagonal(ones(p)),
                            Diagonal(ones(p)),
                            true,
                            false)
end

function make_xml(ml::MBDXMLElement)
    ml.el = make_MBD(ml.precision, ml.diagonal_prec)

    pm_el = ml.el[bn.PRECISION_MATRIX][1]

    if ml.is_random
        mp_el = pm_el[bn.MATRIX_PARAMETER][1]

        ml.prior_el = make_Wishart_prior(ml.prior_scale, mp_el,
                                        bn.DEFAULT_MBD_PRIOR)
    end
end

function make_MBD(prec::AbstractArray{Float64, 2},
        diagonal_prec::Bool)

    dim = size(prec, 1)
    el = new_element(bn.MULTIVARIATE_DIFFUSION_MODEL)
    set_attribute(el, bn.ID, bn.DIFFUSION_ID)
    prec_el = new_element(bn.PRECISION_MATRIX)
    if diagonal_prec
        @assert isdiag(prec)
        add_diagonal_matrix(prec_el, diag(prec), id = bn.DIFF_PREC_ID, lower="0")
    else
        add_matrix_parameter(prec_el, prec, id = bn.DIFF_PREC_ID)
    end
    add_child(el, prec_el)
    return el
end

function get_loggables(el::MBDXMLElement)
    if !el.is_random
        return Vector{XMLElement}(undef, 0)
    end
    make_xml(el)
    return [el.el[bn.PRECISION_MATRIX][1][bn.MATRIX_PARAMETER][1]]

end

function get_matrix_parameter(el::MBDXMLElement)
    return get_loggable(el)
end

function get_priors(xml::MBDXMLElement)
    if !xml.is_random
        return Vector{XMLElement}(undef, 0)
    end
    make_xml(xml)
    return [xml.prior_el]
end
