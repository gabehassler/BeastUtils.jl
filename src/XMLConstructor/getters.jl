# getting elements


function get_data(bx::BEASTXMLElement)
    return find_element(bx, DataXMLElement)
end

function get_newick(bx::BEASTXMLElement)
    return find_element(bx, NewickXMLElement)
end

function get_treeModel(bx::BEASTXMLElement)
    return find_element(bx, TreeModelXMLElement)
end

function get_mbd(bx::BEASTXMLElement)
    return find_element(bx, MBDXMLElement)
end

function get_extension(bx::BEASTXMLElement)
    return find_element(bx, ModelExtensionXMLElement)
end

function get_repeatedMeasures(bx::BEASTXMLElement)
    return find_element(bx, RepeatedMeasuresXMLElement)
end

function get_integratedFactorModel(bx::BEASTXMLElement)
    return find_element(bx, IntegratedFactorsXMLElement)
end

function get_traitLikelihood(bx::BEASTXMLElement)
    return find_element(bx, TraitLikelihoodXMLElement)
end

function get_operators(bx::BEASTXMLElement)
    ops = find_element(bx, OperatorsXMLElement)
    return ops.els
end

function get_mcmc(bx::BEASTXMLElement)
    return find_element(bx, MCMCXMLElement)
end

function get_timer(bx::BEASTXMLElement)
    return find_element(bx, TimerXMLElement)
end

function get_loadings_op(bx::BEASTXMLElement)
    ops = get_operators(bx)
    for op in ops
        t = typeof(op)
        if t <: LoadingsGibbsOperatorXMLElement
            return op
        elseif t <: HMCOperatorXMLElement
            for grad in op.grads
                g = typeof(grad)
                if g <: LoadingsGradientXMLElement || g <: FactorLoadingsGradientXMLElement
                    return op
                end
            end
        end
    end
    error("No loadings operator.")
end

function get_multiplicative_gamma_op(bx::BEASTXMLElement)
    ops = get_operators(bx)
    for op in ops
        t = typeof(op)
        if t <: NormalGammaPrecisionOperatorXMLElement
            if typeof(op.ggp) <: MultipilcativeGammaGibbsProvider
                return op
            end
        end
    end
    error("No multiplicative gamma gibbs operator.")
end
