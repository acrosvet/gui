using Pkg

Pkg.activate(".")
Pkg.instantiate()
Pkg.precompile()

Pkg.add("PackageCompiler")

using PackageCompiler

create_sysimage([:CSV, :DataFrames, :Dates, :Distributions, :Genie, :GenieFramework, :GraphPlot, :Graphs, :JLD2, :LatinHypercubeSampling, :LightGraphs, :Logging, :Parameters, :PlotlyBase, :Random, :SparseArrays, :StatsBase], sysimage_path = "installed_libs.so")

exit()