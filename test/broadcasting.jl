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
    @test st == sin.(t)
    @test st.(pi/3) == evaluate(st, pi/3)
    @test st(pi/3) == evaluate.(st, pi/3)
    @test st.([0.0, pi/3]) == evaluate(st, [0.0, pi/3])

    # Nested Taylor1 tests
    t = Taylor1(Int, 3)
    ts = similar(t);
    @test get_order(t) == get_order(ts) == 3
    @test typeof(ts) == Taylor1{Int}
    ts .= t
    @test ts == t
    @. ts = 3 * t^2 - 1
    @test ts == 3 * t^2 - 1

    tt = Taylor1([zero(t), one(t)], 2)
    tts = similar(tt);
    @test typeof(tts) == Taylor1{Taylor1{Int}}
    @test get_order(tt) == get_order(tts) == 2
    @test tt .== tt
    @. tts = 3 * tt^2 - 1
    @test tts == 3 * tt^2 - 1
    @test eltype(similar(1.0*tt)) == Taylor1{Float64}

    ttt = Taylor1([zero(tt), one(tt)])
    ttts = similar(ttt);
    @test typeof(ttts) == Taylor1{Taylor1{Taylor1{Int}}}
    @test get_order(ttt) == get_order(ttts) == 1
    @test ttt .≈ ttt
    @. ttts = 3 * ttt^2 - 1
    @test ttts == -1
end

@testset "Broadcasting for HomogeneousPolynomial and TaylorN" begin
    x, y = set_variables("x y", order=3)
    xH = x[1]
    yH = y[1]

    @test xH .== xH
    @test yH .≈ yH
    @test xH .== xH
    @test x[2] .== y[2]

    xHs = similar(xH)
    @test typeof(xHs) == typeof(xH)
    @test get_order(xHs) == get_order(xH)
    xHs .= xH
    @test xHs == xH
    @. xHs = 2 * xH + yH
    @test xHs == 2 * xH + yH

    @test 1 .* xH == xH
    @test 1 .* [xH] == [xH]
    @test [1] .* xH == [xH]

    @test x .== x
    @test y .≈ y
    @test x .!= (1 + x)

    p = similar(x)
    @test typeof(p) == typeof(x)
    @test get_order(p) == get_order(x)
    p .= x
    @test p == x
    @. p = 1 + 2*x + 3x^2 - x * y
    @test p == 1 + 2*x + 3*x^2 - x * y

    @test 1.0 .+ x == 1.0 + x
    @test y .+ x == y + x
    @test typeof(big"1.0" .+ x) == TaylorN{BigFloat}
    @test 1.0 .+ [x] == [1.0 + x]
    @test typeof(1.0 .+ [y]) == Vector{TaylorN{Float64}}
    @test 1.0 .+ [x, 2y] == [1.0 + x, 1.0 + 2y]
    @test [1.0,2.0] .+ [x 2y] == [1.0+x 1.0+2y; 2.0+x 2.0+2y]
    @test [1.0] .+ x == x .+ [1.0] == [1.0 + x]
    @test 1.0 .* y == y
    @test typeof(1.0 .* x .* y) == TaylorN{Float64}
end
