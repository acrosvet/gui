# Generate a network graph

using LightGraphs
using Graphs
using GraphPlot

function make_static_graph(nfarms, nedges, in_degree, out_degree )
    g = Graphs.SimpleGraphs.static_scale_free(nfarms, nedges, in_degree, out_degree, rng = MersenneTwister(42), finite_size_correction = true)

    space = DataFrame(Graphs.edges(g))

    @save "data/static_network.jld2" space

   # return space
end

