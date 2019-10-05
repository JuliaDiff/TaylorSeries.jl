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
    "fateman40.jl"
    )

for file in testfiles
    # Skipping tests with intervals, since IntervalArithmetic.jl requires Julia v1.1+
    VERSION < v"1.1" && file == "intervals.jl" && continue
    include(file)
end
