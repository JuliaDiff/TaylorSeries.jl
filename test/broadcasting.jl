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
    ts.coeffs .= t.coeffs
    @test ts == t
    ts.coeffs .= getcoeff.(3 * t^2 - 1, 0:order(t))
    @test ts == 3 * t^2 - 1

    #  `tt` has to be `Taylor1{Taylor1{Float64}}` (instead of `Taylor1{Taylor1{Int}}`)
    # since the method a^n (n integer) is equivalent to `a^float(n).`
    tt = Taylor1([zero(1.0*t), one(t)], 2)
    tts = zero(tt)
    @test tt .== tt
    tts.coeffs .= getcoeff.(3 * tt^2 - 1, 0:order(tt))
    @test tts == 3 * tt^2 - 1

    ttt = Taylor1([zero(tt), one(tt)])
    ttts = zero(ttt)
    @test ttt .≈ ttt
    ttts.coeffs .= getcoeff.(3 * ttt^1 - 1, 0:order(ttt))
    @test ttts == 3 * ttt^1 - 1
    ttts.coeffs .= getcoeff.(3 * ttt^3 - 1, eachindex(ttts))
    @test ttts == - 1.0
end

@testset "Broadcasting with HomogeneousPolynomial and TaylorN" begin
    x, y = variables!("x y", order=3, nowarn=true)
    xH = x[1]
    yH = y[1]

    @test xH .== xH
    @test yH .≈ yH
    @test xH .== xH
    @test x[2] .== y[2]

    xHs = zero(xH)
    xHs.coeffs .= xH.coeffs
    @test xHs == xH
    xHs.coeffs .= getfield(2 * xH + yH, :coeffs)
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
    p.coeffs .= x.coeffs
    @test p == x
    p.coeffs .= getfield(1 + 2*x + 3x^2 - x * y, :coeffs)
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

@testset "Broadcasting with mixtures Taylor1{TaylorN{T}}" begin
    x, y = variables!("x", numvars=2, order=6, nowarn=true)
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
    tNs.coeffs .= tN.coeffs
    @test tNs == tN
    tNs.coeffs .= getfield(y[1] * tN^2 - 1, :coeffs)
    @test tNs == y[1] * tN^2 - 1
    tNs.coeffs .= getfield(y * tN^2 - 1, :coeffs)
    @test tNs == y * tN^2 - 1
end

@testset "Broadcasting with mixtures TaylorN{Taylor1{T}}" begin
    variables!("x", numvars=2, order=6, nowarn=true)
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

    # Note: Broadcasting assignement changed in Julia v1.14; see #444
    # Use `copyto!` with nested Taylor1's
    tN1s = zero(tN1)
    # tN1s .= tN1
    copyto!(tN1s.coeffs, tN1.coeffs)
    @test tN1s == tN1
    # @. tN1s = t * tN1^2 - 1
    copyto!(tN1s.coeffs, getproperty(t * tN1^2 - 1, :coeffs))
    @test tN1s == t * tN1^2 - 1
end
