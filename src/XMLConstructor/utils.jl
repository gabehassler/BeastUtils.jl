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
    mbd_el = bx.MBD_el
    rm_el = bx.extension_el

    diffVar_el = MatrixInverseXMLElement(mbd_el)
    rmVar_el = MatrixInverseXMLElement(rm_el)

    diffCor_el = CorrelationMatrixXMLElement(mbd_el, true)
    rmCor_el = CorrelationMatrixXMLElement(rm_el, true)

    vp_el = VarianceProportionXMLElement(bx.traitLikelihood_el, bx.treeModel_el,
                        bx.extension_el, bx.MBD_el)

    loggables = LoggablesXMLElement([diffVar_el, rmVar_el, diffCor_el, rmCor_el, vp_el],
                                [false, false, false, false, false])

    if isnothing(bx.loggables)
        bx.loggables = loggables
    else
        bx.loggables = join_loggables(bx.loggables, loggables)
    end
    if !isnothing(bx.mcmc_el)
        bx.mcmc_el.loggables = bx.loggables
    end

    return loggables
end

function use_dates!(bx::BEASTXMLElement)
    bx.data_el.use_times = true
    bx.newick_el.fix_tree = false
end
