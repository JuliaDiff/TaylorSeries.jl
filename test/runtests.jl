# This file is part of TaylorSeries.jl, MIT licensed
#
# Tests for TaylorSeries

testfiles = ("onevariable.jl", "manyvariables.jl", "mixtures.jl",
"identities_Euler.jl", "fateman40.jl")

@async for file in testfiles
    include(file)
end
