using BeastUtils.RTrees

newick = "(((A:0.5,B:0.2):0.1,C:0.3):1.0,(D:0.3,E:0.4):0.7);"
tree = RTrees.parse_newick(newick)

taxa1 = ["A", "B", "C", "D", "E"]

Ψ1 = RTrees.tree_variance(tree, taxa1)

Ψ1_true = [
    1.6 1.1 1.0 0.0 0.0;
    1.1 1.3 1.0 0.0 0.0;
    1.0 1.0 1.3 0.0 0.0;
    0.0 0.0 0.0 1.0 0.7;
    0.0 0.0 0.0 0.7 1.1
    ]

tol = 1e-14

diff1 = maximum(abs.(Ψ1 - Ψ1_true))

@test diff1 < tol

perm = [5, 2, 3, 1, 4]
taxa2 = taxa1[perm]
Ψ2_true = Ψ1_true[perm, perm]
Ψ2 = RTrees.tree_variance(tree, taxa2)

diff2 = maximum(abs.(Ψ2 - Ψ2_true))

@test diff2 < tol
