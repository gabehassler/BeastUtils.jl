function df_to_matrix(df::DataFrame) #Converts a data frame with missing values to a matrix with NaN
    taxa = Vector{String}(df[!, :taxon])
    n, p = size(df)
    p = p - 1
    nms = names(df)

    data = fill(NaN, n, p)
    for i = 1:p
        for j = 1:n
            x = df[j, nms[i + 1]]
            if !ismissing(x)
                data[j, i] = x
            end
        end
    end
    return taxa, data
end

function add_MBD_loggables!(bx::BEASTXMLElement)
    mbd_el = get_mbd(bx)
    rm_el = get_repeatedMeasures(bx)
    like_el = get_traitLikelihood(bx)
    treeModel_el = get_treeModel(bx)

    diffVar_el = MatrixInverseXMLElement(mbd_el)
    rmVar_el = MatrixInverseXMLElement(rm_el)

    diffCor_el = CorrelationMatrixXMLElement(mbd_el, true)
    rmCor_el = CorrelationMatrixXMLElement(rm_el, true)

    vp_el = VarianceProportionXMLElement(like_el, treeModel_el, rm_el, mbd_el)

    loggables = LoggablesXMLElement([diffVar_el, rmVar_el, diffCor_el, rmCor_el, vp_el],
                                [false, false, false, false, false])

    add_loggables(bx, loggables)

    return loggables
end

function use_dates!(bx::BEASTXMLElement)
    data = get_data(bx)
    newick = get_newick(bx)
    data.use_times = true
    newick.fix_tree = false
end

function add_trait!(bx::BEASTXMLElement, data::Matrix{Float64}, trait::String)
    p = size(data, 2)
    data_el = get_data(bx)
    tm_el = get_treeModel(bx)

    add_trait!(data_el, data, trait)
    add_leaf_param!(tm_el, trait, p)
end
