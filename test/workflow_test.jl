using Test
using BeastUtils
using BeastUtils.TreeUtils, BeastUtils.Simulation, BeastUtils.MatrixUtils,
      BeastUtils.SlowLikelihoods, BeastUtils.XMLConstructor, BeastUtils.Logs,
      BeastUtils.RunBeast
using PhyloNetworks
import Random
Random.seed!(666)

const tol = 1e-10

N = 10
K = 2
P = 5

taxa = ["taxon_" * "$i" for i = 1:N]

tree = rtree(N, labels = taxa)
newick = writeTopology(tree)

## Testing Basic MBD Model
Σ = randpd(P)

tdm = TreeDiffusionModel(tree, Σ)
tsm = TraitSimulationModel(taxa, tdm)

data = simulate(tsm)


ll = loglikelihood(tsm, data)

# TODO: move above to separate file


## Testing MBD Models
Σ = randpd(P)
Γ = randpd(P)

tdm = TreeDiffusionModel(tree, Σ)
rvm = ResidualVarianceModel(Γ)
tsm = TraitSimulationModel(taxa, tdm, rvm)

data = simulate(tsm)


bx = make_residual_xml(data, taxa, newick)
set_mbd_precision(bx, inv(Σ))
set_residual_precision(bx, inv(Γ))


filename = "test_residual"
xml_path = "$filename.xml"
save_xml(xml_path, bx)
run_beast(xml_path)
log_path = "$filename.log"
col, likelihoods = get_log_match(log_path, "likelihood", burnin=0.0)

standardize_data!(data)
ll_slow = loglikelihood(tsm, data, rescale_tree = true)
ll_beast = likelihoods[1]

@test ll_slow ≈ ll_beast atol=tol