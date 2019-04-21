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
    @test [1.0] .+ t == 1.0 + t  # This should return a vector !!
    @test t .+ t == 2t == 2 .* t

    st = sin(t)
    @test st(pi/3) == evaluate(st, pi/3)
    @test st.(pi/3) == evaluate(st, pi/3)
    @test st.([0.0, pi/3]) == evaluate(st, [0.0, pi/3])

end

@testset "Broadcasting for nested Taylor1's" begin
    t = Taylor1(Int, 3)
    @test typeof(similar(t)) == Taylor1{Int}
    @test get_order(t) == 3

    tt = Taylor1([zero(t), one(t), zero(t)])
    @test typeof(similar(tt)) == Taylor1{Taylor1{Int}}
    @test get_order(tt) == get_order(similar(tt)) == 2

    ttt = Taylor1([zero(tt), one(tt)])
    @test typeof(similar(ttt)) == Taylor1{Taylor1{Taylor1{Int}}}
    @test get_order(ttt) == get_order(similar(ttt)) == 1

    tf = Taylor1([zero(1.0*t), one(t)])
    @test typeof(similar(tf)) == Taylor1{Taylor1{Float64}}
    @test get_order(tf) == get_order(similar(tf)) == 1
end
