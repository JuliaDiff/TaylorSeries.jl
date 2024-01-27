# This file is part of TaylorSeries.jl, MIT licensed
#
# Tests for TaylorSeries

testfiles = (
    "aqua.jl",
    "onevariable.jl",
    "manyvariables.jl",
    "mixtures.jl",
    "mutatingfuncts.jl",
    "intervals.jl",
    "broadcasting.jl",
    "identities_Euler.jl",
    "fateman40.jl"
    )

for file in testfiles
    include(file)
end
# After `using intervalArithmetic` new ambiguities may arise
include("aqua.jl")
