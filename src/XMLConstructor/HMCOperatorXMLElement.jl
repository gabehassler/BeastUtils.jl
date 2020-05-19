abstract type AbstractGradientXMLElement <: MyXMLElement end


mutable struct HMCOperatorXMLElement <: OperatorXMLElement
    el::XMLOrNothing
    param_provider::MyXMLElement
    grads::Array{AbstractGradientXMLElement}
    weight::Real
    n_steps::Int
    step_size::Float64
    draw_variance::Float64
    auto_optimize::Bool
    check_grad::Int
    grad_tolerance::Float64

    function HMCOperatorXMLElement(param_provider::MyXMLElement,
                        grads::Array{AbstractGradientXMLElement})
        return new(nothing, param_provider, grads,
                    1.0, 10, 0.05, 1.0, true, 0, 1e-3)
    end

end

function make_xml(hmcxml::HMCOperatorXMLElement)
    el = new_element(bn.HMC)
    set_attributes(el, [bn.WEIGHT => string(hmcxml.weight),
                        bn.NSTEPS => string(hmcxml.n_steps),
                        bn.STEP_SIZE => string(hmcxml.step_size),
                        bn.DRAW_VARIANCE => string(hmcxml.draw_variance),
                        bn.AUTO_OPTIMIZE => string(hmcxml.auto_optimize),
                        bn.GRAD_CHECK_COUNT => string(hmcxml.check_grad),
                        bn.GRAD_TOLERANCE => string(hmcxml.grad_tolerance)])

    jg_el = new_child(el, bn.JOINT_GRADIENT)

    for grad in hmcxml.grads
        make_xml(grad)
        add_child(jg_el, grad.el)
    end

    add_ref_el(el, get_hmc_parameter(hmcxml.param_provider))

    hmcxml.el = el
    return el
end

mutable struct LoadingsGradientXMLElement <: AbstractGradientXMLElement
    el::XMLOrNothing
    ifxml::IntegratedFactorsXMLElement

    function LoadingsGradientXMLElement(ifxml::IntegratedFactorsXMLElement)

        return new(nothing, ifxml)
    end

end

function make_xml(gxml::LoadingsGradientXMLElement)
    ifxml = gxml.ifxml
    make_xml(ifxml)

    if isnothing(ifxml.msls)

        el = new_element(bn.GRADIENT)
        add_ref_el(el, get_normal_prior(gxml.ifxml))
        add_ref_el(el, gxml.ifxml.loadings_el)


    else
        make_xml(ifxml.msls, ifxml.loadings_el)
        el = reference_element(ifxml.msls.ms_el)
    end

    gxml.el = el
    return el
end

mutable struct FactorLoadingsGradientXMLElement <: AbstractGradientXMLElement
    el::XMLOrNothing
    ifxml::IntegratedFactorsXMLElement
    tdlxml::TraitLikelihoodXMLElement

    function FactorLoadingsGradientXMLElement(
            ifxml::IntegratedFactorsXMLElement,
            tdlxml::TraitLikelihoodXMLElement)

        return new(nothing, ifxml, tdlxml)
    end

end

function make_xml(flgxml::FactorLoadingsGradientXMLElement)
    make_xml(flgxml.ifxml)
    make_xml(flgxml.tdlxml)

    el = new_element(bn.INTEGRATED_FACTOR_GRADIENT)
    add_ref_el(el, flgxml.ifxml.el)
    add_ref_el(el, flgxml.tdlxml.el)

    flgxml.el = el
    return el
end










# <hamiltonianMonteCarloOperator weight="5.0" nSteps="1" stepSize="0.05" drawVariance="1.0" autoOptimize="true" gradientCheckCount="100000">
#             <jointGradient>
#                 <!--                <gradient>-->
#                 <!--                    <distributionLikelihood idref="L.prior"/>-->
#                 <!--                    <matrixParameter idref="L"/>-->
#                 <!--                </gradient>-->
#                 <loadingsShrinkageGradient id="L.grad">
#                     <matrixParameter idref="L"/>
#                     <rowPriors>
#                         <bayesianBridge idref="bb1"/>
#                         <bayesianBridge idref="bb2"/>
#                     </rowPriors>
#                 </loadingsShrinkageGradient>
#                 <integratedFactorAnalysisLoadingsGradient threadType="parallel">
#                     <integratedFactorModel idref="factorModel"/>
#                     <traitDataLikelihood idref="traitLikelihood"/>
#                 </integratedFactorAnalysisLoadingsGradient>
#             </jointGradient>
#             <matrixParameter idref="L"/>
#         </hamiltonianMonteCarloOperator>
