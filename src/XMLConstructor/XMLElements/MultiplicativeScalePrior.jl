mutable struct MultiplicativeScalePrior <: MyXMLElement
    mults::Parameter
    precs::MultiplicativeParameter
    mult_prior::MultivariateGammaLikelihood
    scale_prior::NormalMatrixNormLikelihood
end

function make_xml(mps::MultiplicativeScalePrior)
    make_xml(mps.mults)
    make_xml(mps.precs)
    make_xml(mps.mult_prior)
    make_xml(mps.scale_prior)
    return nothing
end

function get_loggables(mps::MultiplicativeScalePrior)
    make_xml(mps)
    return [mps.mults.el]
end

function NormalGammaPrecisionOperatorXMLElement(
    multiplicativeLikelihood::MultiplicativeScalePrior)


    ggp = MultipilcativeGammaGibbsProvider(multiplicativeLikelihood.mults,
                multiplicativeLikelihood.scale_prior)

    return NormalGammaPrecisionOperatorXMLElement(nothing,
                ggp,
                multiplicativeLikelihood.mult_prior,
                1.0)
end