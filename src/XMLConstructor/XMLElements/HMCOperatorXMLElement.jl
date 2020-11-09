mutable struct HMCOperatorXMLElement <: OperatorXMLElement
    el::XMLOrNothing
    param_provider::MyXMLElement
    grads::Array{<:MyXMLElement}
    already_made::Vector{Bool}
    weight::Real
    n_steps::Int
    step_size::Float64
    draw_variance::Float64
    auto_optimize::Bool
    check_grad::Int
    grad_tolerance::Float64
    geodesic::Bool
    transform::String
    mask::Vector{<:Real}

    function HMCOperatorXMLElement(param_provider::MyXMLElement,
                        grads::Array{<:MyXMLElement};
                        geodesic::Bool = false,
                        already_made::Vector{Bool} = fill(false, length(grads)),
                        transform::String = "")
        return new(nothing, param_provider, grads, already_made,
                    1.0, 5, 0.05, 1.0, true, 0, 1e-3, geodesic, transform,
                    Float64[])
    end

end

function make_xml(hmcxml::HMCOperatorXMLElement)
    name = hmcxml.geodesic ? bn.GEODESIC_HMC : bn.HMC
    el = new_element(name)
    set_attributes(el, [bn.WEIGHT => string(hmcxml.weight),
                        bn.NSTEPS => string(hmcxml.n_steps),
                        bn.STEP_SIZE => string(hmcxml.step_size),
                        bn.DRAW_VARIANCE => string(hmcxml.draw_variance),
                        bn.AUTO_OPTIMIZE => string(hmcxml.auto_optimize),
                        bn.GRAD_CHECK_COUNT => string(hmcxml.check_grad),
                        bn.GRAD_TOLERANCE => string(hmcxml.grad_tolerance)])

    jg_el = new_child(el, bn.JOINT_GRADIENT)

    for i = 1:length(hmcxml.grads)
        grad = hmcxml.grads[i]
        make_xml(grad)
        if hmcxml.already_made[i]
            add_ref_el(jg_el, grad.el)
        else
            add_child(jg_el, grad.el)
        end
    end

    add_ref_el(el, get_hmc_parameter(hmcxml.param_provider))
    if hmcxml.transform != ""
        t_el = new_child(el, bn.TRANSFORM)
        set_attribute(t_el, bn.TYPE, hmcxml.transform)
    end

    if length(hmcxml.mask) > 0
        m_el = new_child(el, bn.MASK)
        add_child(m_el, make_parameter(value=hmcxml.mask))
    end


    hmcxml.el = el
    return el
end

mutable struct LoadingsGradientXMLElement <: MyXMLElement
    el::XMLOrNothing
    ifxml::IntegratedFactorsXMLElement

    function LoadingsGradientXMLElement(ifxml::IntegratedFactorsXMLElement)

        return new(nothing, ifxml)
    end

end

function make_xml(gxml::LoadingsGradientXMLElement)
    ifxml = gxml.ifxml
    el = make_loadings_gradient(ifxml.loadings, ifxml.loadings_prior)

    gxml.el = el
    return el
end

mutable struct FactorLoadingsGradientXMLElement <: MyXMLElement
    el::XMLOrNothing
    ifxml::IntegratedFactorsXMLElement
    tdlxml::TraitLikelihoodXMLElement
    name::String


    function FactorLoadingsGradientXMLElement(
            ifxml::IntegratedFactorsXMLElement,
            tdlxml::TraitLikelihoodXMLElement)

        return new(nothing, ifxml, tdlxml, bn.INTEGRATED_FACTOR_GRADIENT)
    end

end

function NormalizedLoadingsGradientXMLElement(
    ifxml::IntegratedFactorsXMLElement,
    tdlxml::TraitLikelihoodXMLElement)
    grad = FactorLoadingsGradientXMLElement(ifxml, tdlxml)
    grad.name = bn.NORMALIZED_LOADINGS_GRADIENT
    return grad
end

function ScaleLoadingsGradientXMLElement(
    ifxml::IntegratedFactorsXMLElement,
    tdlxml::TraitLikelihoodXMLElement)
    grad = FactorLoadingsGradientXMLElement(ifxml, tdlxml)
    grad.name = bn.SCALE_LOADINGS_GRADIENT
    return grad
end

function make_xml(flgxml::FactorLoadingsGradientXMLElement)
    make_xml(flgxml.ifxml)
    make_xml(flgxml.tdlxml)

    el = new_element(flgxml.name)
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
