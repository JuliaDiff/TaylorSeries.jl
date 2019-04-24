# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries

using Test

@testset "Broadcasting for Taylor1 expansions" begin
    t = Taylor1(Int, 5)

    # @test t .= t
    @test t .== t
    @test t .≈ t
    @test t .!= (1 + t)

    @test 1.0 .+ t == 1.0 + t
    @test typeof(1.0 .+ t) == Taylor1{Float64}
    @test 1.0 .+ [t] == [1.0 + t]
    @test typeof(1.0 .+ [t]) == Vector{Taylor1{Float64}}
    @test 1.0 .+ [t, 2t] == [1.0 + t, 1.0 + 2t]
    @test [1.0,2.0] .+ [t 2t] == [1.0+t 1.0+2t; 2.0+t 2.0+2t]
    @test [1.0] .+ t == [1.0 + t]
    @test 1.0 .* t == t
    @test typeof(1.0 .* t) == Taylor1{Float64}

    st = sin(t)
    @test st .== st
    @test st(pi/3) == evaluate(st, pi/3)
    @test st.(pi/3) == evaluate(st, pi/3)
    @test st.([0.0, pi/3]) == evaluate(st, [0.0, pi/3])
end

@testset "Broadcasting for HomogeneousPolynomial and TaylorN" begin
    x, y = set_variables("x y", order=3)
    xH = x[1]
    yH = y[1]

    @test xH .== xH
    @test yH .≈ yH
    @test x[2] .== y[2]

    @test x .== x
    @test y .≈ y
    @test x .!= (1 + x)

    @test 1.0 .+ x == 1.0 + x
    @test y .+ x == y + x
    @test typeof(big"1.0" .+ x) == TaylorN{BigFloat}
    @test 1.0 .+ [x] == [1.0 + x]
    @test typeof(1.0 .+ [y]) == Vector{TaylorN{Float64}}
    @test 1.0 .+ [x, 2y] == [1.0 + x, 1.0 + 2y]
    @test [1.0,2.0] .+ [x 2y] == [1.0+x 1.0+2y; 2.0+x 2.0+2y]
    @test [1.0] .+ x == [1.0 + x]
    @test 1.0 .* y == y
    @test typeof(1.0 .* x .* y) == TaylorN{Float64}
    #
    # @test 1.0 .+ t == 1.0 + t
    # @test typeof(1.0 .+ t) == Taylor1{Float64}
    # @test 1.0 .+ [t] == [1.0 + t]
    # @test typeof(1.0 .+ [t]) == Vector{Taylor1{Float64}}
    # @test 1.0 .+ [t, 2t] == [1.0 + t, 1.0 + 2t]
    # @test [1.0,2.0] .+ [t 2t] == [1.0+t 1.0+2t; 2.0+t 2.0+2t]
    # @test [1.0] .+ t == [1.0 + t]
    # @test 1.0 .* t == t
    # @test typeof(1.0 .* t) == Taylor1{Float64}

    # xH = HomogeneousPolynomial([one(t), zero(t)])
    # @assert typeof(similar(xHt1)) == HomogeneousPolynomial{Taylor1{Int64}}
end

@testset "Broadcasting for nested Taylor1's" begin
    t = Taylor1(Int, 3)
    # @test typeof(similar(t)) == Taylor1{Int}
    @test get_order(t) == 3

    tt = Taylor1([zero(t), one(t), zero(t)])
    # @test typeof(similar(tt)) == Taylor1{Taylor1{Int}}
    # @test get_order(tt) == get_order(similar(tt)) == 2
    @test tt .== tt

    ttt = Taylor1([zero(tt), one(tt)])
    # @test typeof(similar(ttt)) == Taylor1{Taylor1{Taylor1{Int}}}
    # @test get_order(ttt) == get_order(similar(ttt)) == 1
    @test ttt .≈ ttt

    tf = Taylor1([zero(1.0*t), one(t)])
    # @test typeof(similar(tf)) == Taylor1{Taylor1{Float64}}
    # @test get_order(tf) == get_order(similar(tf)) == 1
    @test tf .== tt
end
