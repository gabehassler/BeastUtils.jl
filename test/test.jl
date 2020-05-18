using BeastUtils.TreeUtils
using PhyloNetworks

using Random
# Random.seed!(666)

x = "((A:1.0,B:1.0)E:1.0,C:1.0);"
x2 = "(A:1.0,B:1.0,C:1.0);"

net = readTopology(x2)

t = TreeUtils.rtree(5, ultrametric = true, keep_order = true)
n = writeTopology(t)
println(n)
println(string.(names(vcv(t))))
println()
clipboard("n <- \"$n\"")

# plot(t, showEdgeLength=true)


pbm = ParamsBM(0.0, 1.0)
x = simulate(t, pbm)
