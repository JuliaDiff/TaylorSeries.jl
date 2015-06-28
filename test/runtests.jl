# Tests for TaylorSeries implementation
using TaylorSeries
using FactCheck
using Compat

FactCheck.setstyle(:compact)
# FactCheck.onlystats(true)

include("Taylor1_tests.jl")
include("TaylorN_tests.jl")
include("identities_tests.jl")



exitstatus()
