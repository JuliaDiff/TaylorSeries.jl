# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries

using Test

@testset "Broadcasting for Taylor1 expansions" begin
    t = Taylor1(Int, 5)

    @test 1.0 .+ t == 1.0 + t
    @test typeof(1.0 .+ t) == Taylor1{Float64}
    @test 1.0 .+ [t] == [1.0 + t]
    @test typeof(1.0 .+ [t]) == Vector{Taylor1{Float64}}
    @test 1.0 .+ [t, 2t] == [1.0 + t, 1.0 + 2t]
    @test [1.0,2.0] .+ [t 2t] == [1.0+t 1.0+2t; 2.0+t 2.0+2t]
    @test [1.0] .+ t == 1.0 + t
end
