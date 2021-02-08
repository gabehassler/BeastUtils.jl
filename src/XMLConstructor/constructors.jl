# contstructors for specific, commonly run analyses

function make_pfa_xml(data::Matrix{Float64}, taxa::Vector{T},
            newick::String, k::Int;
            chain_length::Int=100,
            useHMC::Bool=true,
            timing::Bool=false,
            log_factors::Bool=false,
            # shrink_loadings::Bool=false,
            fle::Int=10,
            sle::Int=100) where T <: AbstractString

    beastXML = BEASTXMLElement()
    data_el = DataXMLElement(data, taxa, newick)
    add_child(beastXML, data_el)

    newick_el = NewickXMLElement(newick)
    add_child(beastXML, newick_el)

    treeModel_el = TreeModelXMLElement(newick_el, size(data, 2))
    add_child(beastXML, treeModel_el)

    mbd_el = MBDXMLElement(k, diagonal_prec=true)
    mbd_el.is_random = false
    add_child(beastXML, mbd_el)

    if_el = IntegratedFactorsXMLElement(treeModel_el, k)
    add_child(beastXML, if_el)

    # if shrink_loadings
    #     if_el.loadings_prior = MatrixShrinkageLikelihoods(
    #                     get_loadings_param(
    #                         get_integratedFactorModel(beastXML)
    #                     )
    #                     )
    # end

    traitLikelihood_el = TraitLikelihoodXMLElement(mbd_el, treeModel_el, if_el)
    add_child(beastXML, traitLikelihood_el)

    traitLikelihood_el.attrs[bn.TRAIT_NAME] = bn.DEFAULT_FACTOR_NAME
    traitLikelihood_el.attrs[bn.ALLOW_SINGULAR] = bn.TRUE
    traitLikelihood_el.attrs[bn.STANDARDIZE] = bn.FALSE

    loadings_op = LoadingsGibbsOperatorXMLElement(if_el,
                                                    traitLikelihood_el)
    if useHMC
        prior_grad = LoadingsGradientXMLElement(if_el)
        like_grad = FactorLoadingsGradientXMLElement(if_el,
                                                        traitLikelihood_el)

        loadings_op = HMCOperatorXMLElement(if_el, [prior_grad, like_grad])
    end

    normal_gamma_op = NormalGammaPrecisionOperatorXMLElement(if_el,
                                                                traitLikelihood_el)

    ops_vec = [loadings_op, normal_gamma_op]
    # if shrink_loadings
    #     push!(ops_vec, ShrinkageScaleOperators(if_el.loadings_prior, if_el))
    # end

    operators_el = OperatorsXMLElement(ops_vec)
    add_child(beastXML, operators_el)

    if log_factors
        traitLog_el = TraitLoggerXMLElement(treeModel_el, traitLikelihood_el)
        add_child(beastXML, traitLog_el)
    end

    mcmc_el = MCMCXMLElement(traitLikelihood_el,
                                mbd_el,
                                if_el,
                                operators_el,
                                chain_length=chain_length)
    add_child(beastXML, mcmc_el)

    mcmc_el.file_logEvery = fle
    mcmc_el.screen_logEvery = sle

    if log_factors
        add_loggable(mcmc_el.loggables, traitLog_el)
    end

    if timing

        mcmc_el.attrs[bn.FULL_EVALUATION] = "0"
        filename = mcmc_el.filename
        timer_el = TimerXMLElement(mcmc_el)
        add_child(beastXML, timer_el)
        timer_el.filename = "$(filename)_timer.txt"
    end

    return beastXML
end

function make_orthogonal_pfa_xml(data::Matrix{Float64}, taxa::Vector{T},
                                 newick::String, k::Int;
                                 shrinkage::Float64 = 1e1,
                                 shrink_first::Bool = false,
                                 chain_length::Int=100,
                                 timing::Bool=false,
                                 log_factors::Bool=false,
                                 rotate_prior::Bool=false,
                                 force_ordered::Bool=false,
                                 fle::Int=10,
                                 sle::Int=100) where T <: AbstractString

    beastXML = BEASTXMLElement()
    data_el = DataXMLElement(data, taxa, newick)
    add_child(beastXML, data_el)

    newick_el = NewickXMLElement(newick)
    add_child(beastXML, newick_el)

    treeModel_el = TreeModelXMLElement(newick_el, size(data, 2))
    add_child(beastXML, treeModel_el)

    mbd_el = MBDXMLElement(k, diagonal_prec=true)
    mbd_el.is_random = false
    add_child(beastXML, mbd_el)

    if_el = IntegratedFactorsXMLElement(treeModel_el, k, orthonormal = true)
    add_child(beastXML, if_el)

    if force_ordered
        load_scale = get_loadings_scale(if_el)
        fol_el = ForceOrderedLikelihood(load_scale)
        add_child(beastXML, fol_el)
    end

    mult_vals = [1.0 for i = 1:k]
    mult_shapes = [shrinkage for i = 1:k]
    mult_scales = ones(k)

    if !shrink_first
        mult_vals[1] = 1.0
        mult_shapes[1] = 1.0
    end

    mults = Parameter(mult_vals, "mults")

    mults_prior = MultivariateGammaLikelihood(mults, mult_shapes, mult_scales,
                                              "mults.likelihood")



    precs = MultiplicativeParameter(mults, "globalPrecision")

    loadings_prior = NormalMatrixNormLikelihood(precs, if_el.loadings, "scale.prior")

    multiplicative_prior = MultiplicativeScalePrior(mults, precs, mults_prior,
                                loadings_prior)

    if_el.loadings_prior = multiplicative_prior

    traitLikelihood_el = TraitLikelihoodXMLElement(mbd_el, treeModel_el, if_el)
    add_child(beastXML, traitLikelihood_el)

    traitLikelihood_el.attrs[bn.TRAIT_NAME] = bn.DEFAULT_FACTOR_NAME
    traitLikelihood_el.attrs[bn.ALLOW_SINGULAR] = bn.TRUE
    traitLikelihood_el.attrs[bn.STANDARDIZE] = bn.FALSE

    U_grad = NormalizedLoadingsGradientXMLElement(if_el, traitLikelihood_el)
    loadings_op = HMCOperatorXMLElement(if_el, [U_grad], geodesic=true)

    scale_grad = ScaleLoadingsGradientXMLElement(if_el, traitLikelihood_el)
    scale_op = HMCOperatorXMLElement(if_el.loadings.scale,
                                     [scale_grad, loadings_prior],
                                     already_made = [false, true],
                                     transform = "log")

    mults_op = NormalGammaPrecisionOperatorXMLElement(if_el.loadings_prior)

    normal_gamma_op = NormalGammaPrecisionOperatorXMLElement(if_el,
                                                        traitLikelihood_el)

    # scale_op = ScaleOperator(if_el.loadings.scale)
    # scale_op.weight = 5.0


    ops_vec = [loadings_op, scale_op, mults_op, normal_gamma_op]

    if force_ordered
        cso = ColumnSwapOperator(if_el.loadings.U, 1.0)
        push!(ops_vec, cso)
    end


    operators_el = OperatorsXMLElement(ops_vec)
    add_child(beastXML, operators_el)

    if log_factors
        traitLog_el = TraitLoggerXMLElement(treeModel_el, traitLikelihood_el)
        add_child(beastXML, traitLog_el)
    end

    mcmc_el = MCMCXMLElement(traitLikelihood_el,
                        mbd_el,
                        if_el,
                        operators_el,
                        chain_length=chain_length)
    if force_ordered
        add_prior!(mcmc_el, fol_el)
    end

    add_child(beastXML, mcmc_el)

    mcmc_el.file_logEvery = fle
    mcmc_el.screen_logEvery = sle

    if log_factors
        add_loggable(mcmc_el.loggables, traitLog_el)
    end

    if timing

        mcmc_el.attrs[bn.FULL_EVALUATION] = "0"
        filename = mcmc_el.filename
        timer_el = TimerXMLElement(mcmc_el)
        add_child(beastXML, timer_el)
        timer_el.filename = "$(filename)_timer.txt"
    end

    return beastXML
end

function make_residual_xml(data::Matrix{Float64}, taxa::Vector{T},
            newick::String; chain_length=100) where T <: AbstractString

    beastXML = BEASTXMLElement()
    data_el = DataXMLElement(data, taxa, newick)
    add_child(beastXML, data_el)

    newick_el = NewickXMLElement(newick)
    add_child(beastXML, newick_el)

    treeModel_el = TreeModelXMLElement(newick_el, size(data, 2))
    add_child(beastXML, treeModel_el)

    MBD_el = MBDXMLElement(size(data, 2))
    add_child(beastXML, MBD_el)

    extension_el = RepeatedMeasuresXMLElement(treeModel_el)
    add_child(beastXML, extension_el)


    traitLikelihood_el = TraitLikelihoodXMLElement(MBD_el, treeModel_el,
                                                    extension_el)
    add_child(beastXML, traitLikelihood_el)


    compoundPrecisionOperator =
        CompoundPrecisionOperatorXMLElement(traitLikelihood_el, MBD_el,
                                            extension_el)

    operators_el = OperatorsXMLElement(compoundPrecisionOperator)
    add_child(beastXML, operators_el)

    mcmc_el = MCMCXMLElement(traitLikelihood_el,
                                MBD_el,
                                extension_el,
                                operators_el,
                                chain_length=chain_length
                                )
    add_child(beastXML, mcmc_el)

    timer_el = TimerXMLElement(mcmc_el)
    add_child(beastXML, timer_el)

    return beastXML
end


function make_joint_xml(newick::String, dm::DataModel, jpm::JointProcessModel)
    beastXML = BEASTXMLElement()
    data_el = DataXMLElement(dm, newick)
    add_child(beastXML, data_el)

    newick_el = NewickXMLElement(newick)
    add_child(beastXML, newick_el)

    treeModel_el = TreeModelXMLElement(newick_el, dm.trait_names,
                                       trait_dimensions(dm))
    add_child(beastXML, treeModel_el)

    mbd_el = MBDXMLElement(jpm.diffusion_model)
    mbd_el.precision = LKJPrecisionXMLElement(tip_dimension(jpm))
    mbd_el.precision_prior = LKJPrecisionPriors(mbd_el.precision)
    add_child(beastXML, mbd_el)

    for i = 1:length(jpm.extensions)
        ext_el = make_xmlelement(jpm.extensions[i], treeModel_el, ind=i)
        add_child(beastXML, ext_el)
    end

    # traitLikelihood_el = TraitLikelihoodXMLElement(MBD_el, treeModel_el,
    #                                                 extension_el)

    return beastXML
end

function make_old_pfa_xml(data::Matrix{Float64}, taxa::Vector{T},
                        newick::String, k::Int;
                        chain_length::Int = 100,
                        timing::Bool = false) where T <: AbstractString

    beastXML = BEASTXMLElement()

    data_el = DataXMLElement([data, zeros(size(data, 1), k)],
                [bn.DEFAULT_TRAIT_NAME, bn.FACTOR_TRAIT_NAME], taxa, newick)
    add_child(beastXML, data_el)

    newick_el = NewickXMLElement(newick)
    add_child(beastXML, newick_el)

    treeModel_el = TreeModelXMLElement(newick_el, size(data, 2))
    add_child(beastXML, treeModel_el)

    MBD_el = MBDXMLElement(k, diagonal_prec=true)
    MBD_el.is_random = false
    add_child(beastXML, MBD_el)

    traitLikelihood_el = TraitLikelihoodXMLElement(MBD_el,
                                                   treeModel_el, nothing)
    traitLikelihood_el.xml_name = bn.MULTIVARIATE_TRAIT_LIKELIHOOD
    traitLikelihood_el.attrs[bn.TRAIT_NAME] = bn.DEFAULT_FACTOR_NAME
    delete!(traitLikelihood_el.attrs, bn.ALLOW_SINGULAR)
    add_child(beastXML, traitLikelihood_el)

    lfm_el = LatentFactorModelXMLElement(treeModel_el,
                                         traitLikelihood_el, k)
    add_child(beastXML, lfm_el)

    loadings_op = OldLoadingsGibbsOperatorXMLElement(lfm_el)
    prec_op = LatentFactorModelPrecisionOperatorXMLElement(lfm_el)
    fac_op = FactorTreeGibbsOperatorXMLElement(lfm_el,
                                               traitLikelihood_el)


    operators_el = OperatorsXMLElement([loadings_op, prec_op, fac_op])
    add_child(beastXML, operators_el)


    mcmc_el = MCMCXMLElement(traitLikelihood_el,
                             MBD_el,
                             lfm_el,
                             operators_el,
                             chain_length = chain_length)
    add_child(beastXML, mcmc_el)

    if timing
        set_full_eval!(mcmc_el, 0)
        filename = mcmc_el.filename
        timer_el = TimerXMLElement(mcmc_el)
        timer_el.filename = "$(filename)_timer.txt"
        add_child(beastXML, timer_el)
    end

    return beastXML
end


################################################################################
## Lower level constructors
################################################################################

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

# ## TODO: need to update below


# function make_oldPFA_XML(data::Matrix{Float64}, taxa::Vector{T},
#     newick::String, k::Int;
#     chain_length::Int = 100,
#     timing::Bool = false) where T <: AbstractString

#     beastXML = BEASTXMLElement()
#     beastXML.data_el = DataXMLElement([data, zeros(size(data, 1), k)],
#                 [bn.DEFAULT_TRAIT_NAME, bn.FACTOR_TRAIT_NAME], taxa, newick)
#     beastXML.newick_el = NewickXMLElement(newick)
#     beastXML.treeModel_el = TreeModelXMLElement(beastXML.newick_el,
#                                                 size(data, 2))
#     beastXML.MBD_el = MBDXMLElement(k)
#     beastXML.MBD_el.is_random = false
#     beastXML.MBD_el.diagonal_prec = true

#     beastXML.traitLikelihood_el = TraitLikelihoodXMLElement(beastXML.MBD_el,
#                         beastXML.treeModel_el, nothing)

#     beastXML.traitLikelihood_el.xml_name = bn.MULTIVARIATE_TRAIT_LIKELIHOOD

#     beastXML.extension_el = LatentFactorModelXMLElement(beastXML.treeModel_el,
#                                             beastXML.traitLikelihood_el, k)


#     beastXML.traitLikelihood_el.attrs[bn.TRAIT_NAME] = bn.DEFAULT_FACTOR_NAME

#     delete!(beastXML.traitLikelihood_el.attrs, bn.ALLOW_SINGULAR)





#     loadings_op = OldLoadingsGibbsOperatorXMLElement(beastXML.extension_el)
#     prec_op = LatentFactorModelPrecisionOperatorXMLElement(beastXML.extension_el)
#     fac_op = FactorTreeGibbsOperatorXMLElement(beastXML.extension_el,
#                 beastXML.traitLikelihood_el)


#     beastXML.operators_el = OperatorsXMLElement([loadings_op, prec_op, fac_op])


#     beastXML.mcmc_el = MCMCXMLElement(beastXML.traitLikelihood_el,
#                                     beastXML.MBD_el,
#                                     beastXML.extension_el,
#                                     beastXML.operators_el,
#                                     chain_length = chain_length)

#     if timing

#     beastXML.mcmc_el.attrs[bn.FULL_EVALUATION] = "0"
#     filename = beastXML.mcmc_el.filename
#     beastXML.timer_el = TimerXMLElement(beastXML.mcmc_el)
#     beastXML.timer_el.filename = "$(filename)_timer.txt"
#     end

#     return beastXML


# end

# function make_validation_MBD_XML(mis_data::Matrix{Float64},
#                                 obs_data::Matrix{Float64},
#                                 taxa::Vector{T},
#                                 newick::String;
#                                 chain_length = 100) where T <: AbstractString

#     @assert size(obs_data) == size(mis_data)

#     beastXML = BEASTXMLElement()
#     beastXML.data_el = DataXMLElement([mis_data, obs_data],
#             [bn.DEFAULT_TRAIT_NAME, "$(bn.DEFAULT_TRAIT_NAME)True"],
#             taxa,
#             newick)

#     beastXML.newick_el = NewickXMLElement(newick)

#     p = size(mis_data, 2)
#     beastXML.treeModel_el = TreeModelXMLElement(
#             beastXML.newick_el,
#             beastXML.data_el.trait_names,
#             [p, p],
#             [bn.LEAF_TRAITS, "$(bn.LEAF_TRAITS)True"]
#             )

#     beastXML.MBD_el = MBDXMLElement(p)
#     beastXML.extension_el = RepeatedMeasuresXMLElement(beastXML.treeModel_el, beastXML.MBD_el)
#     beastXML.extension_el.standardize_traits = false
#     beastXML.traitLikelihood_el = TraitLikelihoodXMLElement(beastXML.MBD_el, beastXML.treeModel_el, beastXML.extension_el)
#     compoundPrecisionOperator = CompoundPrecisionOperatorXMLElement(beastXML.traitLikelihood_el, beastXML.MBD_el, beastXML.extension_el)
#     beastXML.operators_el = OperatorsXMLElement(compoundPrecisionOperator)

#     # traitValidation_el = TraitValidationXMLElement(
#     #             beastXML.treeModel_el,
#     #             beastXML.traitLikelihood_el
#     #             )
#     #
#     # traitValidation_el.standardize = beastXML.extension_el.standardize_traits
#     #
#     # validation_el = CrossValidationXMLElement(traitValidation_el)
#     model_extension_el = ModelExtensionLoggerXMLElement(beastXML.traitLikelihood_el,
#                                                 beastXML.extension_el)

#     beastXML.loggables = LoggablesXMLElement(
#             [model_extension_el],
#             [false]
#             )

#     add_MBD_loggables!(beastXML)

#     beastXML.mcmc_el = MCMCXMLElement(beastXML.traitLikelihood_el,
#             beastXML.MBD_el,
#             beastXML.extension_el,
#             beastXML.operators_el,
#             beastXML.loggables,
#             chain_length = chain_length)

#     beastXML.timer_el = TimerXMLElement(beastXML.mcmc_el)

#     return beastXML
#     end

#     function make_timing_MBD_XML(data::Matrix{Float64}, taxa::Vector{String},
#                         newick::String; chain_length::Int = 1)

#     bx = make_MBD_XML(data, taxa, newick, chain_length = chain_length)
#     fpo_el = FireParameterOperatorXMLElement(bx.extension_el)
#     bx.operators_el = OperatorsXMLElement(fpo_el)
#     bx.mcmc_el.operators = bx.operators_el
#     bx.mcmc_el.log_files = false
#     bx.mcmc_el.attrs[bn.FULL_EVALUATION] = "0"
#     filename = bx.mcmc_el.filename
#     bx.timer_el.filename = "$(filename)_timer.txt"
#     return bx
# end


# function make_xml_oldPFA(bx::BEASTXMLElement)
#     xdoc = XMLDocument()
#     bx.el = new_element(bn.BEAST)
#     set_root(xdoc, bx.el)
#     add_child(bx.el, make_xml(bx.data_el))
#     add_child(bx.el, make_xml(bx.newick_el))
#     add_child(bx.el, make_xml(bx.treeModel_el))
#     add_el(bx, bx.MBD_el)
#     add_el(bx, bx.traitLikelihood_el)
#     add_el(bx, bx.extension_el)
#     add_el(bx, bx.operators_el)
#     add_el(bx, bx.varProp_el)
#     add_el(bx, bx.misc_els)
#     add_loggables(bx)
#     add_el(bx, bx.mcmc_el)
#     add_el(bx, bx.timer_el)
#     return xdoc
# end
