module XMLTools

using LightXML, Statistics, LinearAlgebra

export collect_data,
    data_stats,
    remove_nans,
    replace_data!,
    introduce_sparsity!,
    get_corrs,
    complete_cases!,
    replace_newick!,
    replace_filelog!,
    replace_chainLength!,
    replace_logEvery!,
    remove_all_data!,
    remove_dates!,
    get_element_by_id,
    standardize_data!,
    replace_all_data!,
    get_trait_dimmension,
    add_note!,
    transfer_data!,
    get_newick,
    set_diffusion_precision!,
    set_rm_precision!,
    set_operator_weight!,
    setup_trait_validation!,
    multivariateDistributionLikelihood,
    leafTraitParameter,
    insert_element,
    get_likelihood,
    get_diffusion_precision,
    get_rm_precision,
    get_extension_precision,
    get_loadings_mat,
    get_if_precision,
    collect_dates,
    make_matrix_parameter

const TAXA = "taxa"
const TAXON = "taxon"
const ID = "id"
const NAME = "name"
const MISSING_VAL = "NA"
const NEWICK = "newick"
const MCMC = "mcmc"
const LOG = "log"
const FILENAME = "fileName"
const CHAINLENGTH = "chainLength"
const LOGEVERY = "logEvery"
const DATE = "date"
const OPERATORS = "operators"
const ATTR = "attr"
const MULTIVARIATE_DIFFUSION_MODEL = "multivariateDiffusionModel"
const PRECISION_MATRIX = "precisionMatrix"
const MATRIX_PARAMETER = "matrixParameter"
const PARAMETER = "parameter"
const VALUE = "value"
const SAMPLING_PRECISION = "samplingPrecision"
const REPEATED_MEASURES = "repeatedMeasuresModel"
const WEIGHT = "weight"
const BEAST = "beast"
const TRAIT_VALIDATION = "traitValidation"
const TRAIT_NAME = "traitName"
const INFERRED_TRAIT = "inferredTrait"
const USE_TREE_TRAITS = "useTreeTraits"
const TRAIT_DATA_LIKELIHOOD = "traitDataLikelihood"
const IDREF = "idref"
const TRAIT_PARAMETER = "traitParameter"
const FILELOG = "fileLog"
const MULT_DIST_LIKELIHOOD = "multivariateDistributionLikelihood"
const DISTRIBUTION = "distribution"
const MVM_DISTRIBUTION = "multivariateNormalDistributionModel"
const MEAN_PARAMETER = "meanParameter"
const PRECISION_PARAMETER = "precisionParameter"
const DATA = "data"
const LEAF_TRAIT_PARAMETER = "leafTraitParameter"
const TREE_MODEL = "treeModel"
const LIKELIHOOD = "likelihood"
const POSTERIOR = "posterior"
const INTEGRATED_FACTOR_MODEL = "integratedFactorModel"
const PRECISION = "precision"
const LOADINGS = "loadings"
const DIAGONAL_MATRIX = "DiagonalMatrix"
const DIMENSION = "dimension"
const DATE = "date"

function get_child_by_attribute(xel::XMLElement, attr::String, attr_value::String)
    children = child_elements(xel)
    for child in children
        if attribute(child, attr) == attr_value
            return child
        end
    end
    throw("ERROR: No children have attribute $attr with value \"$attr_value\".")
end

function get_element_by_id(el::XMLElement, id::String)
    children = child_elements(el)
    final = nothing
    for child in children
        if attribute(child, ID) == id
            final = child
            break
        else
            potential = get_element_by_id(child, id)
            if potential != nothing
                final = potential
                break
            end
        end
    end
    return final
end

function get_element_by_id(xdoc::XMLDocument, id::String)
    el = root(xdoc)
    return get_element_by_id(el, id)
end


function collect_data(xdoc::XMLDocument, trait_name::String)
    taxa_el = find_element(root(xdoc), TAXA)
    taxa_children = collect(child_elements(taxa_el))
    n = length(taxa_children)
    taxa = Vector{String}(undef, n)
    p = length(split(content(get_child_by_attribute(taxa_children[1], NAME,
            trait_name))))
    data = Matrix{Float64}(undef, n, p)
    for i = 1:n
        taxa[i] = attribute(taxa_children[i], ID)
        taxa_vals = split(content(get_child_by_attribute(taxa_children[i],
                NAME, trait_name)))
        for j = 1:p
            if taxa_vals[j] == MISSING_VAL
                data[i, j] = NaN
            else
                data[i, j] = parse(Float64, taxa_vals[j])
            end
        end
    end

    return taxa, data
end

function collect_data(path::String, trait_name::String)
    xdoc = parse_file(path)
    to_return = collect_data(xdoc, trait_name)
    free(xdoc)
    return to_return
end

function collect_dates(xdoc::XMLDocument)
    xroot = root(xdoc)
    taxa_el = find_element(xroot, TAXA)
    taxa_children = collect(child_elements(taxa_el))
    n = length(taxa_children)
    dates = zeros(n)
    taxa = Vector{String}(undef, n)
    for i = 1:n
        taxa[i] = attribute(taxa_children[i], ID)

        date_el = find_element(taxa_children[i], DATE)
        dates[i] = parse(Float64, attribute(date_el, VALUE))
    end
    return taxa, dates
end


function collect_dates(path::String)
    xdoc = parse_file(path)
    to_return = collect_dates(xdoc)
    free(xdoc)
    return to_return
end

function replace_data!(xdoc::XMLDocument, trait_name::String,
        replacement_data::Dict{String, Vector{Float64}})
    taxa_el = find_element(root(xdoc), TAXA)
    for taxon_el in child_elements(taxa_el)
        taxon = attribute(taxon_el, ID, required=true)
        trait_el = get_child_by_attribute(taxon_el, NAME, trait_name)
        set_content(trait_el, make_data_string(replacement_data[taxon]))
    end
end

function remove_all_data!(xdoc::XMLDocument)
    taxa_el = find_element(root(xdoc), TAXA)
    unlink(taxa_el)
    free(taxa_el)
end

function transfer_data!(src_doc::XMLDocument, dest_doc::XMLDocument)
    src_taxa = find_element(root(src_doc), TAXA)
    dest_taxa = find_element(root(dest_doc), TAXA)
    for node in child_nodes(dest_taxa)
        unlink(node)
        free(node)
    end
    for node in child_nodes(src_taxa)
        add_child(dest_taxa, node)
    end
end


function remove_dates!(xdoc::XMLDocument)
    taxa_el = find_element(root(xdoc), TAXA)
    for taxon in child_elements(taxa_el)
        date = find_element(taxon, DATE)
        unlink(date)
        free(date)
    end
end


function introduce_sparsity!(xdoc::XMLDocument, trait_name::String,
        sparsity::Vector{Float64})

    child_els = child_elements(find_element(root(xdoc), TAXA))
    n  = length(collect(child_els))
    int_sparsity = [Int(floor(x * n)) for x in sparsity]
    return introduce_sparsity!(xdoc, trait_name, int_sparsity)
end

function introduce_sparsity!(xdoc::XMLDocument, trait_name::String,
        sparsity::Vector{Int})

    taxa_el = find_element(root(xdoc), TAXA)
    taxa = collect(child_elements(taxa_el))
    n = length(taxa)
    removed_inds = get_removed_inds(sparsity, n)

    p = length(sparsity)
    removed_taxa = [Vector{String}(undef, sparsity[i]) for i = 1:p]
    for j = 1:p
        ind = 0
        for i in removed_inds[j]
            ind += 1
            taxon_el = taxa[i]
            trait_el = get_child_by_attribute(taxon_el, NAME, trait_name)
            data = split(content(trait_el))
            data[j] = MISSING_VAL
            set_content(trait_el, join(data, ' '))
            removed_taxa[j][ind] = attribute(taxon_el, ID, required = true)
        end
    end
    return removed_taxa
end

function get_removed_inds(num_removed::Vector{Int}, n::Int)
    p = length(num_removed)
    for i = 1:p
        @assert num_removed[i] <= n
    end
    removed_inds = [Vector{Int}(undef, num_removed[i]) for i = 1:p]
    for i = 1:p
        for j = 1:num_removed[i]
            ok = false
            while !ok
                draw = rand(1:n)
                if !(draw in removed_inds[i])
                    removed_inds[i][j] = draw
                    ok = true
                end
            end
        end
        sort!(removed_inds[i])
    end
    return removed_inds
end




function make_data_string(data::Vector{Float64}; missing_val = MISSING_VAL)
    p = length(data)
    s_data = Vector{String}(undef, p)
    for i = 1:p
        if isnan(data[i])
            s_data[i] = missing_val
        else
            s_data[i] = string(data[i])
        end
    end
    return join(s_data, " ")
end


function data_stats(data::Matrix{Float64})
    n, p = size(data)
    n_observed = zeros(Int, p)
    mean_sums = zeros(p)
    var_sums = zeros(p)
    for i = 1:n
        for j = 1:p
            if !isnan(data[i, j])
                n_observed[j] += 1
                mean_sums[j] += data[i, j]
                var_sums[j] += data[i, j]^2
            end
        end
    end

    for j = 1:p
        mean_sums[j] = mean_sums[j] / n_observed[j]
        var_sums[j] = var_sums[j] / n_observed[j] - mean_sums[j]^2
    end
    cormat = get_corrs(data)
    perc_missing = [1 - n_observed[i] / n for i = 1:p]
    return mean_sums, var_sums, cormat, perc_missing
end

function remove_nans(x::Vector{Float64})
    removed = Vector{Float64}()
    for n in x
        if !isnan(n)
            push!(removed, n)
        end
    end
    return removed
end

function remove_nans(x1::Vector{Float64}, x2::Vector{Float64})
    removed = Matrix{Float64}(undef, 0, 2)
    n = length(x1)
    for i = 1:n
        if !isnan(x1[i]) && !isnan(x2[i])
            removed = [removed; [x1[i] x2[i]]]
        end
    end
    return removed
end

function get_corrs(data::Matrix{Float64})
    n, p = size(data)
    cormat = zeros(p, p)
    for i = 1:n
        for j = i:p
            pairwise = remove_nans(data[:, i], data[:, j])
            cormat[i, j] = cor(pairwise[:, 1], pairwise[:, 2])
        end
    end
    return cormat
end

function complete_cases!(xdoc::XMLDocument, trait_name::String)
    taxa_el = find_element(root(xdoc), TAXA)
    removed_taxa = Vector{String}(undef, 0)
    for taxon_el in child_elements(taxa_el)
        trait_el = get_child_by_attribute(taxon_el, NAME, trait_name)
        data = split(content(trait_el))
        if MISSING_VAL in data
            unlink(taxon_el)
            push!(removed_taxa, attribute(taxon_el, ID, required = true))
        end
    end
    return removed_taxa
end

function replace_newick!(xdoc::XMLDocument, new_newick::AbstractString)
    newick_el = find_element(root(xdoc), NEWICK)
    strip(new_newick)
    strip(new_newick, ';')
    set_content(newick_el, "$new_newick;")
end

function get_newick(xdoc::XMLDocument)
    newick_el = find_element(root(xdoc), NEWICK)
    newick = content(newick_el)
    newick = strip(newick)
    # newick = strip(newick, ';')
    return String(newick)
end

function get_newick(path::String)
    xdoc = parse_file(path)
    newick = get_newick(xdoc)
    free(xdoc)
    return newick
end

function replace_filelog!(xdoc::XMLDocument, new_filename::String)
    mcmc_el = find_element(root(xdoc), MCMC)
    changed = false
    for element in child_elements(mcmc_el)
        if name(element) == LOG && has_attribute(element, FILENAME)
            set_attribute(element, FILENAME, new_filename)
            changed = true
            break
        end
    end
    @assert changed
end

function replace_chainLength!(xdoc::XMLDocument, new_chainLength::Int)
    mcmc_el = find_element(root(xdoc), MCMC)
    set_attribute(mcmc_el, CHAINLENGTH, string(new_chainLength))
end

function replace_logEvery!(xdoc::XMLDocument, new_logEvery::Int)
    mcmc_el = find_element(root(xdoc), MCMC)
    changed = false
    for element in child_elements(mcmc_el)
        if name(element) == LOG && has_attribute(element, LOGEVERY) &&
                has_attribute(element, FILENAME)
            set_attribute(element, LOGEVERY, string(new_logEvery))
            changed = true
            break
        end
    end
    @assert changed
end

function standardize_data!(xdoc::XMLDocument, trait::String;
        keep_zeros::Bool = false)
    taxa, data = collect_data(xdoc, trait)
    n, p = size(data)
    standardize_data!(data, keep_zeros = keep_zeros)
    d = Dict{String, Vector{Float64}}()
    for i = 1:n
        d[taxa[i]] = data[i, :]
    end
    replace_data!(xdoc, trait, d)
end

function standardize_data!(data::Matrix{Float64}; keep_zeros::Bool = false)
    n, p = size(data)
    means, vars = missing_MVs(data, keep_zeros = keep_zeros)
    for j = 1:p
        for i = 1:n
            if !(keep_zeros && data[i, j] == 0.0)
                data[i, j] = (data[i, j] - means[j]) / sqrt(vars[j])
            end
        end
    end
end


function missing_means(X::Array{Float64, 2}; keep_zeros::Bool = false)
    n, p = size(X)
    means = Vector{Float64}(undef, p)
    for j = 1:p
        t = 0.0
        N = 0
        for i = 1:n
            if !isnan(X[i, j])
                if !(keep_zeros && X[i, j] == 0.0)
                    t += X[i, j]
                    N += 1
                end
            end
        end
        means[j] = t / N
    end
    return means
end

function missing_vars(X::Matrix{Float64};
        means::Vector{Float64} = Vector{Float64}(undef, 0),
        keep_zeros::Bool = false)

    if length(means) == 0
        means = missing_means(X)
    end
    n, p = size(X)
    vars = Vector{Float64}(undef, p)
    for j = 1:p
        t = 0.0
        N = 0
        for i = 1:n
            if !isnan(X[i, j])
                if !(keep_zeros && X[i, j] == 0.0)
                    t += (X[i, j] - means[j])^2
                    N += 1
                end
            end
        end
        vars[j] = t / (N - 1)
    end
    return vars
end

function missing_MVs(X::Matrix{Float64}; keep_zeros::Bool = false)
    means = missing_means(X, keep_zeros = keep_zeros)
    vars = missing_vars(X, means = means, keep_zeros = keep_zeros)
    return means, vars
end


function replace_all_data!(xdoc::XMLDocument, taxa::Vector{String},
        data::Matrix{Float64};
        traits::Vector{String} = Vector{String}(undef, 0),
        all_name::String = "traits",
        dates::Vector{Float64} = Float64[])

    if length(dates) > 0
        error("not implemented")
    end
    n, p = size(data)
    use_traits = false
    if length(traits) != 0
        use_traits = true
    end
    @assert length(taxa) == n
    if use_traits
        @assert length(traits) == p
    end
    taxa_el = find_element(root(xdoc), TAXA)
    children = child_nodes(taxa_el)
    for child in children
        unlink(child)
        free(child)
    end
    all_traits = zeros(p)
    for i = 1:n

        new_taxon = new_element(TAXON)
        set_attribute(new_taxon, ID, taxa[i])
        if use_traits
            for j = 1:p
                new_attr = new_element(ATTR)
                set_attribute(new_attr, NAME, traits[j])
                set_content(new_attr, string(data[i, j]))
                add_child_level(new_taxon, new_attr, 3)
            end
        end
        for j = 1:p
            all_traits[j] = data[i, j]
        end
        all_attr = new_element(ATTR)
        set_attribute(all_attr, NAME, all_name)
        set_content(all_attr, make_data_string(all_traits))
        add_child_level(new_taxon, all_attr, 3)
        # add_text(taxa_el, "\t\t")
        add_child_level(taxa_el, new_taxon, 2)
        # @show has_children(new_taxon)
        # @show length(collect(child_elements(new_taxon)))
        # add_text(taxa_el, "\n")
    end
end

function replace_all_data!(xdoc::XMLDocument, taxa::Vector{String},
        data::Vector{Array{Float64, 2}}, traits::Vector{String};
        dates::Vector{Float64} = Float64[])

    use_dates = false
    if length(dates) > 0
        @assert length(dates) == length(taxa)
        use_dates = true
    end

    p = length(data)
    @assert length(traits) == p
    n = length(taxa)

    taxa_el = find_element(root(xdoc), TAXA)
    children = child_nodes(taxa_el)
    for child in children
        unlink(child)
        free(child)
    end
    all_traits = zeros(p)
    for i = 1:n
        new_taxon = new_element(TAXON)
        set_attribute(new_taxon, ID, taxa[i])
        if use_dates
            date = new_element(DATE)
            set_attribute(date, VALUE, dates[i])
            set_attribute(date, "direction", "forwards")
            set_attribute(date, "units", "years")
            add_child_level(new_taxon, date, 3)
        end
        for j = 1:p
            new_attr = new_element(ATTR)
            set_attribute(new_attr, NAME, traits[j])
            set_content(new_attr, make_data_string(data[j][i, :]))
            add_child_level(new_taxon, new_attr, 3)
        end
        add_child_level(taxa_el, new_taxon, 2)
    end
end



function get_trait_dimmension(xdoc::XMLDocument, trait::String)
    taxa_el = find_element(root(xdoc), TAXA)
    child = collect(child_elements(taxa_el))[1]
    trait_el = get_child_by_attribute(child, NAME, trait)
    p = length(split(content(trait_el)))
    return p
end

function add_child_level(parent_el::XMLElement, child_el::XMLElement,
        level::Int; line_below::Bool = false)
    add_text(parent_el, "\n")
    tabs = join(["\t" for i = 1:level])
    add_text(parent_el, tabs)
    # @show has_children(child_el)
    if length(collect(child_elements(child_el))) > 0
        add_text(child_el, "\n$tabs")
    end
    add_child(parent_el, child_el)
    if line_below
        add_text(parent_el, "\n")
    end
end

function add_note!(xdoc::XMLDocument, note::String)
    xroot = root(xdoc)
    comment_note = "<!--\n$note\n-->\n"
    t = new_textnode(comment_note)
    add_child(xroot, t)
end

function set_matrix_parameter!(xel::XMLElement, P::Matrix{Float64})
    p = size(P, 1)
    i = 0
    for child in child_elements(xel)
        i += 1
        set_attribute(child, VALUE, join(P[i, :], " "))
    end
    @assert i == p
end

function set_diffusion_precision!(xdoc::XMLDocument,
        P::Matrix{Float64})

    @assert isposdef(P)
    mp_el = get_diffusion_precision_el(xdoc)
    set_matrix_parameter!(mp_el, P)
end

function get_diffusion_precision(xdoc::XMLDocument)
    mp_el = get_diffusion_precision_el(xdoc)
    P = get_matrix_parameter(mp_el)
    return P
end

function get_diffusion_precision_el(xdoc::XMLDocument)
    diffusion_el = find_element(root(xdoc), MULTIVARIATE_DIFFUSION_MODEL)
    pmat_el = find_element(diffusion_el, PRECISION_MATRIX)
    mp_el = find_element(pmat_el, MATRIX_PARAMETER)
    if mp_el == nothing
        mp_el = find_element(pmat_el, DIAGONAL_MATRIX)
    end
    return mp_el
end

function set_rm_precision!(xdoc::XMLDocument, P::Matrix{Float64})

    @assert isposdef(P)
    mp_el = get_rm_precision_el(xdoc)
    set_matrix_parameter!(mp_el, P)
end

function get_rm_precision(xdoc::XMLDocument)
    mp_el = get_rm_precision_el(xdoc)
    return get_matrix_parameter(mp_el)
end

function get_rm_precision_el(xdoc::XMLDocument)
    rm_el = find_element(root(xdoc), REPEATED_MEASURES)
    pmat_el = find_element(rm_el, SAMPLING_PRECISION)
    mp_el = find_element(pmat_el, MATRIX_PARAMETER)
    return mp_el
end

function get_extension_precision(xdoc::XMLDocument)
    xroot = root(xdoc)
    rm_el = find_element(xroot, REPEATED_MEASURES)
    if_el = find_element(xroot, INTEGRATED_FACTOR_MODEL)
    if rm_el == nothing
        if if_el == nothing
            error("Could not find $(REPEATED_MEASURES) or $(INTEGRATED_FACTOR_MODEL) elements under $(name(xroot))")
        else
            return get_if_precision(xdoc) #TODO
        end
    else
        if if_el == nothing
            return get_rm_precision(xdoc)
        else
            error("XML contained both $(REPEATED_MEASURES) and $(INTEGRATED_FACTOR_MODEL) elements. Did not know which one to get precision from.")
        end
    end
end

function get_if_precision(xdoc::XMLDocument)
    if_el = find_element(root(xdoc), INTEGRATED_FACTOR_MODEL)
    prec_el = find_element(if_el, PRECISION)
    param_el = find_element(prec_el, PARAMETER)
    prec_vec = get_vector_parameter(param_el)
    return diagm(0 => prec_vec)
end

function get_loadings_mat(xdoc::XMLDocument)
    xroot = root(xdoc)
    if_el = find_element(xroot, INTEGRATED_FACTOR_MODEL)
    loadings_el = find_element(if_el, LOADINGS)
    mp_el = find_element(loadings_el, MATRIX_PARAMETER)
    if has_children(mp_el)
        return get_matrix_parameter(mp_el)
    else
        id = attribute(mp_el, IDREF)
    end
    mp_el = get_element_by_id(xdoc, id)
    return get_matrix_parameter(mp_el)
end




function set_operator_weight!(xdoc::XMLDocument, operator_name::String,
        weight::Int)
    operator_el = find_element(root(xdoc), OPERATORS)
    xel = find_element(operator_el, operator_name)
    set_attribute(xel, WEIGHT, weight)
end

function insert_element(xdoc::XMLDocument,
        element::XMLElement,
        reference_name::String;
        above::Bool = true)

    # xdoc = parse_file(src_path)

    xdoc2 = XMLDocument()
    create_root(xdoc2, BEAST)
    xroot2 = root(xdoc2)
    original_els = child_elements(root(xdoc))
    n_element = deepcopy(element)
    for el in original_els
        # free(xdoc)
        # free(xdoc2)
        # return 0
        # @show name(el)
        new_el = xml_deepcopy(el)
        if name(el) != reference_name
            add_child_level(xroot2, new_el, 1)
        else
            if above
                add_child_level(xroot2, n_element, 1)
                add_child_level(xroot2, new_el, 1)
            else
                add_child_level(xroot2, new_el, 1)
                add_child_level(xroot2, n_element, 1)
            end
        end

    end
#
    for el in child_nodes(root(xdoc))
        unlink(el)
        free(el)
    end
    for el in child_elements(root(xdoc2))
        add_child(root(xdoc), xml_deepcopy(el))

    end
    free(xdoc2)
    return 0
end

function xml_deepcopy(el::XMLElement)
    xdoc = parse_string(string(el))
    new_el = root(xdoc)
    unlink(new_el)
    free(xdoc)
    return new_el
end

# function setup_trait_validation!(xdoc::XMLDocument;
#         id::String = "validation",
#         true_name::String = "traitsTrue",
#         true_param::String = "leafTraitsTrue",
#         use_mask::Bool = false,
#         use_tree_traits::Bool = false,
#         inferred_trait::String = "")
#
#     xroot = root(xdoc)
#     @show 1
#
#     validation = new_element(TRAIT_VALIDATION)
#     set_attribute(validation, TRAIT_NAME, true_name)
#     set_attribute(validation, "logSum", "false")
#     if use_tree_traits
#         set_attribute(validation, USE_TREE_TRAITS, "true")
#         set_attribute(validation, INFERRED_TRAIT, inferred_trait)
#     end
#     @show 2
#     trait_likelihood_el = pull_element(xdoc, TRAIT_DATA_LIKELIHOOD)
#     add_child(validation, trait_likelihood_el)
#     trait_param_el = new_element(TRAIT_PARAMETER)
#     param_el = new_element(PARAMETER)
#     set_attribute(param_el, IDREF, true_param)
#     add_child_level(trait_param_el, param_el, 3)
#     add_child_level(validation, trait_param_el, 2)
#
#     xdoc = insert_element!(xdoc, validation, MCMC)
#     save_file(xdoc, joinpath(Directories.desktop, "test.xml"))
#
#
#     log_el = get_element_by_id(xdoc, FILELOG)
#     new_el = new_element(TRAIT_VALIDATION)
#     set_attribute(new_el, IDREF, id)
#     add_child(log_el, new_el)
# end

function pull_element(xdoc::XMLDocument, el_name::String;
        parent::XMLElement = root(xdoc))

    el = find_element(parent, el_name)
    id = attribute(el, ID)
    new_el = new_element(el_name)
    set_attribute(new_el, IDREF, id)
    return new_el
end

function multivariateDistributionLikelihood(id::String, mean_param_id::String,
        prec_param_id::String,
        data::Vector{T}) where T <: Real

    el = new_element(MULT_DIST_LIKELIHOOD)
    set_attribute(el, ID, id)
    dist_el = new_element(DISTRIBUTION)
    mvn_el = new_element(MVM_DISTRIBUTION)
    mean_param = new_element(MEAN_PARAMETER)
    param = new_element(PARAMETER)
    set_attribute(param, IDREF, mean_param_id)
    add_child_level(mean_param, param, 4)
    add_child_level(mvn_el, mean_param, 3)
    prec_param = new_element(PRECISION_PARAMETER)
    mat_param = new_element(MATRIX_PARAMETER)
    set_attribute(mat_param, IDREF, prec_param_id)
    add_child_level(prec_param, mat_param, 4)
    add_child_level(mvn_el, prec_param, 3)
    add_child_level(dist_el, mvn_el, 2)
    add_child_level(el, dist_el, 1)

    data_el = new_element(DATA)
    param_el = new_element(PARAMETER)
    set_attribute(param_el, VALUE, join(data, " "))
    add_child_level(data_el, param_el, 2)
    add_child_level(el, data_el, 1)
    return el
end

function leafTraitParameter(id::String, taxon::String, tree::String,
        param::String)

    el = new_element(LEAF_TRAIT_PARAMETER)
    set_attribute(el, ID, id)
    set_attribute(el, TAXON, taxon)
    tree_el = new_element(TREE_MODEL)
    set_attribute(tree_el, IDREF, tree)
    add_child_level(el, tree_el, 2)
    param_el = new_element(PARAMETER)
    set_attribute(param_el, IDREF, param)
    add_child_level(el, param_el, 2)
    return el
end

function get_likelihood(xdoc::XMLDocument)
    mcmc_el = find_element(root(xdoc), MCMC)
    posterior_el = find_element(mcmc_el, POSTERIOR)
    likelihood_el = find_element(posterior_el, LIKELIHOOD)
    return likelihood_el
end

function get_matrix_parameter(mat_el::XMLElement; T::DataType = Float64)
    if name(mat_el) == MATRIX_PARAMETER
        children = collect(child_elements(mat_el))
        n = length(children)
        v1 = get_vector_parameter(children[1], T = T)
        p = length(v1)

        M = Matrix{T}(undef, n, p)
        M[1, :] .= v1
        for i = 2:n
            M[i, :] = get_vector_parameter(children[i], T = T)
        end
        return M
    elseif name(mat_el) == DIAGONAL_MATRIX
        param_el = find_element(mat_el, PARAMETER)
        dim = parse(Int, attribute(param_el, DIMENSION))
        val = parse(Float64, attribute(param_el, VALUE))
        return diagm(0 => fill(val, dim))
    else
        error("Not yet implemented")
    end
end

function get_vector_parameter(vec_el::XMLElement; T::DataType = Float64)
    buffer = Vector{T}(undef, 0)
    return get_vector_parameter(vec_el, buffer)
end

function get_vector_parameter(vec_el::XMLElement,
        buffer::Vector{T}) where T <: Any
    vals_attr = attribute(vec_el, VALUE, required = true)
    vals = split(vals_attr, ' ')
    p = length(vals)
    if length(buffer) == 0
        buffer = Vector{T}(undef, p)
    else
        @assert length(buffer) == p
    end
    for i = 1:p
        buffer[i] = parse(T, vals[i])
    end
    return buffer
end

function make_sub_ids(mat::Matrix{Float64}, id::String)
    n = size(mat, 1)
    return ["$(id).$i" for i = 1:n]
end

function make_matrix_parameter(mat::Matrix{Float64}, id::String;
                            sub_ids::Vector{String} = make_sub_ids(mat, id))
    el = new_element(MATRIX_PARAMETER)
    set_attribute(el, ID, id)
    n = size(mat, 1)
    for i = 1:n
        p_el = new_child(el, PARAMETER)
        set_attribute(p_el, ID, sub_ids[i])
        set_attribute(p_el, VALUE, join(mat[i, :], ' '))
    end
    return string(el)
end











end
