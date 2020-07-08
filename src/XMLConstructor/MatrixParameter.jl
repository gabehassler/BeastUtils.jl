mutable struct Parameter <: MyXMLElement
    el::XMLOrNothing
    val::T where T <: AbstractVector{Float64}
    dim::Int
    lower::Float64
    upper::Float64
    id::String
end

function Parameter(val::Vector{Float64}, id::String)
    return Parameter(nothing, val, length(val), NaN, NaN, id)
end


function make_parameter(;id::String = "",
        value::AbstractArray{T} = Vector{Int8}(undef, 0),
        lower::Float64 = NaN, upper::Float64 = NaN,
        dim::Int = length(value)) where T <: Real

    el = new_element(bn.PARAMETER)
    if id != ""
        set_attribute(el, bn.ID, id)
    end
    if length(value) != 0
        set_attribute(el, bn.VALUE, join(value, ' '))
    end
    if !isnan(lower)
        set_attribute(el, bn.LOWER, lower)
    end
    if !isnan(upper)
        set_attribute(el, bn.UPPER, upper)
    end
    if dim != length(value)
        set_attribute(el, bn.DIMENSION, dim)
    end

    return el
end

function add_parameter(pel::XMLElement; id::String = "",
                       value::AbstractArray{T} = Vector{Int8}(undef, 0),
                       lower::Float64 = NaN, upper::Float64 = NaN,
                       dim::Int = length(value)) where T <: Real

    el = make_parameter(id = id, value = value, lower = lower, upper = upper,
                        dim = dim)
    add_child(pel, el)
    return el
end


function make_xml(param::Parameter)
    el = make_parameter(id = param.id, value = param.val, lower = param.lower,
                        upper = param.upper, dim = param.dim)
    param.el = el
    return el
end


mutable struct MatrixParameter <: MyXMLElement
    el::XMLOrNothing
    mat::Matrix{Float64}
    id::String
    ids::Vector{String}
end

function MatrixParameter(mat::Matrix{Float64}, id::String, ids::Vector{String})
    m, n = size(mat)
    if length(ids) != m
        error("The matrix has $m rows, but $(length(ids)) row ids were provided.")
    end

    return MatrixParameter(nothing, mat, id, ids)
end

function make_xml(mp::MatrixParameter)
    el = new_element(bn.MATRIX_PARAMETER)
    set_attribute(el, bn.ID, mp.id)
    n = size(mp.mat, 1)
    for i = 1:n
        add_parameter(el, id = mp.ids[i], value = @view mp.mat[i, :])
    end
    mp.el = el
    return el
end

function set_mat!(mp::MatrixParameter, mat::S) where S <: AbstractMatrix{Float64}
    if size(mat) != size(mp.mat)
        error("Supplied matrix has dimensions $(size(mat)), while matrix " *
              "parameter has dimensions $(size(mp.mat))")
    end

    mp.mat = mat
end
