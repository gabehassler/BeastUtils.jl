using BeastUtils.XMLConstructor
using BeastUtils.TreeUtils, PhyloNetworks

xc = XMLConstructor

using DataFrames, CSV, Test
### Making sure https://github.com/suchard-group/incomplete_measurements/blob/master/scripts/xml_setup.jl works

function make_mbd(data_path::String, newick_path::String, xml_path::String,
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
        xc.set_data_dates(bx, dates_df[!, :date])
    end
    xc.set_screen_logEvery(bx, 100)
    xc.set_file_logEvery(bx, 10)
    xc.set_filename(bx, filename)
    XMLConstructor.add_MBD_loggables!(bx)



    XMLConstructor.save_xml(xml_path, bx)
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


make_mbd(data_path, newick_path, xml_path, filename)

@test isfile(xml_path)

make_mbd(data_path, newick_path, dates_xml_path, filename, dates_path = dates_path)

@test isfile(dates_xml_path)



### Testing Factor xml
k = 4

# HMC, no shrinkage
bx = XMLConstructor.make_PFA_XML(data, taxa, newick, k, useHMC = false)
XMLConstructor.save_xml("facHMC.xml", bx)
@test isfile("facHMC.xml")

# Gibbs, no shrinkage
bx = XMLConstructor.make_PFA_XML(data, taxa, newick, k, useHMC = true)
XMLConstructor.save_xml("facGibbs.xml", bx)
@test isfile("facGibbs.xml")

# HMC, shrinkage
bx = XMLConstructor.make_PFA_XML(data, taxa, newick, k, useHMC = true,
            shrink_loadings = true)
XMLConstructor.save_xml("facHMCShrink.xml", bx)
@test isfile("facHMCShrink.xml")

# Gibbs, shrinkage
bx = XMLConstructor.make_PFA_XML(data, taxa, newick, k, useHMC = false,
            shrink_loadings = true)
XMLConstructor.save_xml("facGibbsShrink.xml", bx)
@test isfile("facGibbsShrink.xml")
