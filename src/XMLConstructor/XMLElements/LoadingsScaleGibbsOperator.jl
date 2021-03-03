mutable struct LoadingsScaleGibbsOperator <: MyXMLElement
    el::XMLOrNothing
    likelihoods::Vector{MyXMLElement}
    prior::MyXMLElement
    weight::Float64

    function LoadingsScaleGibbsOperator(factors::LatentFactorModelXMLElement,
                                        prior::MyXMLElement)
        return new(nothing, [factors], prior, 1.0)
    end

    function LoadingsScaleGibbsOperator(factors::IntegratedFactorsXMLElement,
                                    tree_likelihood::TraitLikelihoodXMLElement)
        return new(nothing, [factors, tree_likelihood], factors.loadings_prior, 1.0)
    end
end

function make_xml(lsgo::LoadingsScaleGibbsOperator)
    el = new_element(bn.LOCADINGS_SCALE_GIBBS)
    set_attribute(el, bn.WEIGHT, lsgo.weight)
    for like in lsgo.likelihoods
        add_ref_el(el, like)
    end
    add_ref_el(el, get_loadings_prior(lsgo.prior))
    lsgo.el = el
    return el
end
