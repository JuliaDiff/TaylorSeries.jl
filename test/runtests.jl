# This file is part of TaylorSeries.jl, MIT licensed
#
# Tests for TaylorSeries

testfiles = (
    "onevariable.jl",
    "manyvariables.jl",
    "mixtures.jl",
    "mutatingfuncts.jl",
    "identities_Euler.jl",
    "fateman40.jl",
    "intervals.jl"
    )

for file in testfiles
    # Tests in 0.7 pass when current IntervalArithmetics master
    # is used, i.e., IntervalArithmetic 0.14.0, not yet released
    # file == "intervals.jl" && VERSION ≥ v"0.7.0-DEV" && continue
    include(file)
end
