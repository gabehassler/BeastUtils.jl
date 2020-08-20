## set parts of BeastXML

function set_screen_logEvery(bx::BEASTXMLElement, sle::Int)
    mcmc_el = get_mcmc(bx)
    set_screen_logEvery(mcmc_el, sle)
end

function set_file_logEvery(bx::BEASTXMLElement, fle::Int)
    mcmc_el = get_mcmc(bx)
    set_file_logEvery(mcmc_el, fle)
end

function set_filename(bx::BEASTXMLElement, filename::AbstractString)
    mcmc_el = get_mcmc(bx)
    set_filename(mcmc_el, string(filename))
end

function set_data_dates(bx::BEASTXMLElement, dates::AbstractVector{Float64})
    data_el = get_data(bx)
    set_dates(data_el, dates)
end

function set_full_eval(bx::BEASTXMLElement, n_eval::Int)
    mcmc = get_mcmc(bx)
    mcmc.attrs[bn.FULL_EVALUATION] = string(n_eval)
end
