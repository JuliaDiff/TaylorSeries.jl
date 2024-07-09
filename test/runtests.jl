# This file is part of TaylorSeries.jl, MIT licensed
#
# Tests for TaylorSeries

testfiles = (
    "onevariable.jl",
    "manyvariables.jl",
    "mixtures.jl",
    "mutatingfuncts.jl",
    "intervals.jl",
    "broadcasting.jl",
    "identities_Euler.jl",
    "fateman40.jl",
    "staticarrays.jl",
    "jld2.jl",
    "rat.jl",
    # Run Aqua tests at the very end
    "aqua.jl",
    )

for file in testfiles
    include(file)
end
