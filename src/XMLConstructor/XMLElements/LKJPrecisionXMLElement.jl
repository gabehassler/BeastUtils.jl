mutable struct LKJPrecisionXMLElement <: MyXMLElement
    el::XMLOrNothing
    var::AbstractMatrix{Float64}
    diag_param::Parameter
    off_diag_parm::Parameter
    id::String

    function LKJPrecisionXMLElement(var_mat::AbstractArray{Float64},
                                    id::AbstractString)
        check_posdef(var_mat)
        vars, corrs = var_corrs(var_mat)
        diag_param = Parameter(vars, "$id.diagonalVar")
        diag_param.lower = 0.0
        off_diag_param = Parameter(corrs, "$id.correlation")
        off_diag_param.lower = -1.0
        off_diag_param.upper = 1.0
        return new(nothing, var_mat, diag_param, off_diag_param, id)
    end
end

function LKJPrecisionXMLElement(p::Int)
    return LKJPrecisionXMLElement(Diagonal(ones(p)), bn.DIFF_PREC_ID)
end


function make_xml(lkj::LKJPrecisionXMLElement)
    var, corrs = var_corrs(lkj.var)

    set_value(lkj.diag_param, var)
    set_value(lkj.off_diag_parm, corrs)

    el = new_element(bn.CACHED_MATRIX_INVERSE)
    set_id!(el, lkj.id)

    cs_el = new_child(el, bn.COMPOUND_SYMMETRIC_MATRIX)
    set_attributes(cs_el,
                   Dict(bn.AS_CORRELATION => bn.TRUE, bn.IS_CHOLESKY => bn.TRUE)
                  )

    diag_el = new_child(cs_el, bn.DIAGONAL)
    dp_el = make_xml(lkj.diag_param)
    add_child(diag_el, dp_el)
    od_el = new_child(cs_el, bn.OFF_DIAGONAL)
    op_el = make_xml(lkj.off_diag_parm)
    add_child(od_el, op_el)

    lkj.el = el
    return el
end

function var_corrs(var::AbstractArray{Float64})
    vars = diag(var)
    corr_mat = cov2corr(var)
    p = size(var, 1)
    n_corrs = div(p * (p - 1), 2)
    corrs = zeros(n_corrs)
    ind = 0
    for i = 1:p
        for j = (i + 1):p
            ind += 1
            corrs[ind] = corr_mat[i, j]
        end
    end
    return vars, corrs
end


