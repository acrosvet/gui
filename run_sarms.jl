using Pkg

Pkg.activate(".")
Pkg.instantiate()
Pkg.precompile()

using GenieFramework

Genie.loadapp(); Genie.up(async=false)