using BeastUtils.XMLConstructor
using BeastUtils.TreeUtils, PhyloNetworks

using DataFrames, CSV, LightXML, Test
### Making sure https://github.com/suchard-group/incomplete_measurements/blob/master/scripts/xml_setup.jl works

function make_xml(data_path::String, newick_path::String, xml_path::String,
                filename::String; dates_path::String = "")

    df = CSV.read(data_path)

    use_dates = false
    if length(dates_path) > 0
        dates_df = CSV.read(dates_path)
        @assert dates_df[!, :taxon] == df[!, :taxon]
        use_dates = true
    end

    newick = read(newick_path, String)

    taxa, data = XMLConstructor.df_to_matrix(df)


    bx = XMLConstructor.make_MBD_XML(data, taxa, newick, chain_length = 100_000)
    if use_dates
        XMLConstructor.use_dates!(bx)
        bx.data_el.node_times = dates_df[!, :date]
    end
    bx.mcmc_el.screen_logEvery = 100
    bx.mcmc_el.file_logEvery = 10
    bx.mcmc_el.filename = filename
    XMLConstructor.add_MBD_loggables!(bx)


    xdoc = XMLConstructor.make_xml(bx)

    save_file(xdoc, xml_path)
    free(xdoc)
end


n = 10
p = 4
taxa = ["taxon$i" for i = 1:n]
data = randn(n, p)
dates = rand(n)

df = DataFrame()
df.taxon = taxa
for i = 1:p
    df[!, Symbol("trait$i")] = data[:, i]
end

data_path = "data.csv"
dates_path = "dates.csv"
newick_path = "newick.txt"
xml_path = "XMLConstructor.xml"
dates_xml_path = "xml_dates.xml"
filename = "test"


CSV.write(data_path, df)
CSV.write(dates_path, DataFrame(taxon = taxa, date = dates))

newick = writeTopology(TreeUtils.rtree(n, labels = taxa))
write(newick_path, newick)


make_xml(data_path, newick_path, xml_path, filename)

@test isfile(xml_path)

make_xml(data_path, newick_path, dates_xml_path, filename, dates_path = dates_path)

@test isfile(dates_xml_path)
