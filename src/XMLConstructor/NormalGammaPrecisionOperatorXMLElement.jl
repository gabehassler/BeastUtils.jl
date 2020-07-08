abstract type GammaGibbsProvider <: MyXMLElement end

mutable struct NormalGammaPrecisionOperatorXMLElement <: OperatorXMLElement
    el::XMLOrNothing
    ggp::GammaGibbsProvider
    prior_provider::MyXMLElement
    weight::Float64
end

function NormalGammaPrecisionOperatorXMLElement(ggp::GammaGibbsProvider,
                    prior::T) where T <: MyXMLElement
    return NormalGammaPrecisionOperatorXMLElement(nothing, ggp, prior, 1.0)
end


function NormalGammaPrecisionOperatorXMLElement(
            ext_xml::ModelExtensionXMLElement,
            td_xml::TraitLikelihoodXMLElement)
    return NormalGammaPrecisionOperatorXMLElement(nothing,
                NormalExtentsionGibbsProvider(ext_xml, td_xml),
                ext_xml, 1.0)
end

function NormalGammaPrecisionOperatorXMLElement(msl::MatrixShrinkageLikelihoods,
                                                row::Int)
    gpp = MultipilcativeGammaGibbsProvider(msl, row)
    return NormalGammaPrecisionOperatorXMLElement(nothing, gpp, gpp, 1.0)
end


mutable struct NormalExtentsionGibbsProvider <: GammaGibbsProvider
    el::XMLOrNothing
    ext::ModelExtensionXMLElement
    tdl::TraitLikelihoodXMLElement

    function NormalExtentsionGibbsProvider(ext::ModelExtensionXMLElement,
                                           tdl::TraitLikelihoodXMLElement)
        return new(nothing, ext, tdl)
    end
end

function make_xml(negp::NormalExtentsionGibbsProvider)
    el = new_element(bn.NORMAL_EXTENSION)

    make_xml(negp.ext)
    make_xml(negp.tdl)
    set_attribute(el, bn.TREE_TRAIT_NAME, negp.tdl.attrs[bn.TRAIT_NAME])
    add_ref_el(el, negp.ext.el)
    add_ref_el(el, negp.tdl.el)

    negp.el = el
    return el
end

mutable struct MultipilcativeGammaGibbsProvider <: GammaGibbsProvider
    el::XMLOrNothing
    msl::MatrixShrinkageLikelihoods
    row::Int

    function MultipilcativeGammaGibbsProvider(msl::MatrixShrinkageLikelihoods,
                                              row::Int)
        return new(nothing, msl, row)
    end
end

function multiplicative_gamma_gibbs_proivders(msl::MatrixShrinkageLikelihoods)
    k = get_fac_dim(msl)
    mggs = Vector{MultipilcativeGammaGibbsProvider}(undef, k)
    for i = 1:k
        mggs[i] = MultipilcativeGammaGibbsProvider(msl, i)
    end
    return mggs
end

function make_xml(mgg::MultipilcativeGammaGibbsProvider)
    make_xml(mgg.msl)

    el = new_element(bn.MULTIPLICATIVE_GAMMA_GIBBS)
    set_attribute(el, bn.ROW, mgg.row)
    rp_el = new_child(el, bn.ROW_MULTIPLIERS)
    add_ref_el(rp_el, mgg.msl.gp_els[1])
    for mult_el in mgg.msl.mult_els
        add_ref_el(rp_el, mult_el)
    end

    add_ref_el(el, mgg.msl.ms_el)
    mgg.el = el
end

function get_precision_prior(mgg::MultipilcativeGammaGibbsProvider)
    make_xml(mgg.msl)
    return mgg.msl.global_prior_els[mgg.row]
end




function make_xml(ngpxml::NormalGammaPrecisionOperatorXMLElement)
    make_xml(ngpxml.ggp)
    make_xml(ngpxml.prior_provider)

    el = new_element(bn.NGP_OPERATOR)
    set_attribute(el, bn.WEIGHT, ngpxml.weight)
    prior_el = new_child(el, bn.PRIOR)
    add_ref_el(prior_el, get_precision_prior(ngpxml.prior_provider))

    add_child(el, ngpxml.ggp.el)

    ngpxml.el = el
    return el
end
