module gui

using Pkg, CSV, DataFrames, Dates, Distributions, Genie, GenieFramework, GraphPlot, Graphs, JLD2, LatinHypercubeSampling, LightGraphs, Logging, Parameters, PlotlyBase, Random, SparseArrays, StatsBase
#include("../app.jl")
#include("../ui.jl")

Base.@ccallable function julia_main()::Cint
    try
        real_main()
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

function real_main()
   
    Genie.loadapp(); Genie.up(async=false)
end

end