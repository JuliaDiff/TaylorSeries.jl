# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries

using Test

@testset "Broadcasting with Taylor1 expansions" begin
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
    @test [1.0] .+ t == t .+ [1.0] == [1.0 + t]
    @test 1.0 .* t == t
    @test typeof(1.0 .* t) == Taylor1{Float64}

    st = sin(t)
    @test st .== st
    @test st == sin.(t)
    @test st.(pi/3) == evaluate(st, pi/3)
    @test st(pi/3) == evaluate.(st, pi/3)
    @test st.([0.0, pi/3]) == evaluate(st, [0.0, pi/3])

    @test typeof(Float32.(t)) == Taylor1{Float32}
    @test (Float32.(t))[1] == Float32(1.0)
    @test_throws MethodError Float32(t)

    # Nested Taylor1 tests
    t = Taylor1(Int, 3)
    ts = zero(t)
    ts .= t
    @test ts == t
    @. ts = 3 * t^2 - 1
    @test ts == 3 * t^2 - 1

    tt = Taylor1([zero(t), one(t)], 2)
    tts = zero(tt)
    @test tt .== tt
    @. tts = 3 * tt^2 - 1
    @test tts == 3 * tt^2 - 1

    ttt = Taylor1([zero(tt), one(tt)])
    ttts = zero(ttt)
    @test ttt .≈ ttt
    @. ttts = 3 * ttt^2 - 1
    @test ttts == -1
end

@testset "Broadcasting with HomogeneousPolynomial and TaylorN" begin
    x, y = set_variables("x y", order=3)
    xH = x[1]
    yH = y[1]

    @test xH .== xH
    @test yH .≈ yH
    @test xH .== xH
    @test x[2] .== y[2]

    xHs = zero(xH)
    xHs .= xH
    @test xHs == xH
    @. xHs = 2 * xH + yH
    @test xHs == 2 * xH + yH

    @test 1 .* xH == xH
    @test 1 .* [xH] == [xH]
    @test [1] .* xH == xH .* [1] == [xH]

    @test x .== x
    @test y .≈ y
    @test x .!= (1 + x)

    @test typeof(Float32.(x)) == TaylorN{Float32}
    @test (Float32.(x))[1] == HomogeneousPolynomial(Float32[1.0, 0.0])
    @test_throws MethodError Float32(x)

    p = zero(x)
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

@testset "Broadcasting with mixtures Taylor1{TalorN{T}}" begin
    x, y = set_variables("x", numvars=2, order=6)
    tN = Taylor1(TaylorN{Float64}, 3)

    @test tN .== tN
    @test tN .≈ tN
    @test tN .!= (1 + tN)

    @test 1.0 .+ tN == 1.0 + tN
    @test typeof(1.0 .+ tN) == Taylor1{TaylorN{Float64}}
    @test 1.0 .+ [tN] == [1.0 + tN]
    @test typeof(1.0 .+ [tN]) == Vector{Taylor1{TaylorN{Float64}}}
    @test 1.0 .+ [tN, 2tN] == [1.0 + tN, 1.0 + 2tN]
    @test [1.0,2.0] .+ [tN 2tN] == [1.0+tN 1.0+2tN; 2.0+tN 2.0+2tN]
    @test [1.0] .+ tN == tN .+ [1.0] == [1.0 + tN]
    @test 1.0 .* tN == 1.0 * tN
    @test typeof(1.0 .* tN) == Taylor1{TaylorN{Float64}}

    tNs = zero(tN)
    tNs .= tN
    @test tNs == tN
    @. tNs = y[1] * tN^2 - 1
    @test tNs == y[1] * tN^2 - 1
    @. tNs = y * tN^2 - 1
    @test tNs == y * tN^2 - 1
end

@testset "Broadcasting with mixtures TaylorN{Talor1{T}}" begin
    set_variables("x", numvars=2, order=6)
    t = Taylor1(3)
    xHt = HomogeneousPolynomial([one(t), zero(t)])
    yHt = HomogeneousPolynomial([zero(t), t])
    tN1 = TaylorN([HomogeneousPolynomial([t]),xHt,yHt^2])

    @test tN1 .== tN1
    @test tN1 .≈ tN1
    @test tN1 .!= (1 + tN1)

    @test 1.0 .+ tN1 == 1.0 + tN1
    @test typeof(1.0 .+ tN1) == TaylorN{Taylor1{Float64}}
    @test 1.0 .+ [tN1] == [1.0 + tN1]
    @test typeof(1.0 .+ [tN1]) == Vector{TaylorN{Taylor1{Float64}}}
    @test 1.0 .+ [tN1, 2tN1] == [1.0 + tN1, 1.0 + 2tN1]
    @test [1.0, 2.0] .+ [tN1 2tN1] == [1.0+tN1 1.0+2tN1; 2.0+tN1 2.0+2tN1]
    @test [1.0] .+ tN1 == tN1 .+ [1.0] == [1.0 + tN1]
    @test 1.0 .* tN1 == tN1
    @test typeof(1.0 .* tN1) == TaylorN{Taylor1{Float64}}

    tN1s = zero(tN1)
    tN1s .= tN1
    @test tN1s == tN1
    @. tN1s = t * tN1^2 - 1
    @test tN1s == t * tN1^2 - 1
end
