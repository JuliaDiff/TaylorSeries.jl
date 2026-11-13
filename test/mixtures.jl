# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries

using Test
# using LinearAlgebra

@testset "Tests with mixtures of Taylor1 and TaylorN" begin
    set_taylor1_varname(1, "t")
    @test TS.NumberNotSeries == Union{Real,Complex}
    @test TS.NumberNotSeriesN == Union{Real,Complex,Taylor1}

    variables!("x", numvars=2, order=6, nowarn=true)
    xH = HomogeneousPolynomial(Int, 1)
    yH = HomogeneousPolynomial(Int, 2)
    tN = Taylor1(TaylorN{Float64}, 3)
    @test findfirst(tN) == 1
    @test findfirst(zero(tN)) == -1

    @test convert(eltype(tN), tN) == tN
    @test eltype(xH) == HomogeneousPolynomial{Int}
    @test TS.numtype(xH) == Int
    @test eltype(tN) == Taylor1{TaylorN{Float64}}
    @test TS.numtype(tN) == TaylorN{Float64}
    @test normalize_taylor(tN) == tN
    @test order(tN) == 3

    @testset "Lexicographic order: Taylor1{HomogeneousPolynomial{T}} and Taylor1{TaylorN{T}}" begin
        @test HomogeneousPolynomial([1]) > xH > yH > 0.0
        @test -xH^2 < -xH*yH < -yH^2 < -xH^3 < -yH^3 < HomogeneousPolynomial([0.0])
        @test 1 ≥ tN > 2*tN^2 > 100*tN^3 > 0
        @test -2*tN < -tN^2 ≤ 0
    end

    @test string(zero(tN)) == "  0.0 + 𝒪(‖x‖¹) + 𝒪(t⁴)"
    @test string(tN) == " ( 1.0 + 𝒪(‖x‖¹)) t + 𝒪(t⁴)"
    @test string(tN + 3Taylor1(Int, 2)) == " ( 4.0 + 𝒪(‖x‖¹)) t + 𝒪(t³)"
    @test string(xH * tN) == " ( 1.0 x₁ + 𝒪(‖x‖²)) t + 𝒪(t⁴)"

    @test constant_term(xH) == xH
    @test constant_term(tN) == zero(TaylorN([xH]))
    @test linear_polynomial(xH) == xH
    @test linear_polynomial(1+tN+tN^2) == tN
    @test nonlinear_polynomial(1+tN+tN^2) == tN^2

    tN = Taylor1([zero(TaylorN(Float64,1)), one(TaylorN(Float64,1))], 3)
    @test typeof(tN) == Taylor1{TaylorN{Float64}}
    @test string(zero(tN)) == "  0.0 + 𝒪(‖x‖⁷) + 𝒪(t⁴)"
    @test string(tN) == " ( 1.0 + 𝒪(‖x‖⁷)) t + 𝒪(t⁴)"
    @test string(Taylor1([xH+yH])) == "  1 x₁ + 1 x₂ + 𝒪(t¹)"
    @test string(Taylor1([zero(xH), xH*yH])) == " ( 1 x₁ x₂) t + 𝒪(t²)"
    @test string(tN * Taylor1([0,TaylorN([xH+yH])])) == "  0.0 + 𝒪(‖x‖⁷) + 𝒪(t²)"

    t = Taylor1(3)
    xHt = HomogeneousPolynomial(typeof(t), 1)
    @test findfirst(xHt) == 1
    @test convert(eltype(xHt), xHt) === xHt
    @test eltype(xHt) == HomogeneousPolynomial{Taylor1{Float64}}
    @test TS.numtype(xHt) == Taylor1{Float64}
    @test normalize_taylor(xHt) == xHt
    @test string(xHt) == " ( 1.0 + 𝒪(t¹)) x₁"
    xHt = HomogeneousPolynomial([one(t), zero(t)])
    yHt = HomogeneousPolynomial([zero(t), t])
    @test findfirst(yHt) == 2
    @test string(xHt) == " ( 1.0 + 𝒪(t⁴)) x₁"
    @test string(yHt) == " ( 1.0 t + 𝒪(t⁴)) x₂"
    @test string(HomogeneousPolynomial([t])) == " ( 1.0 t + 𝒪(t⁴))"
    @test 3*xHt == HomogeneousPolynomial([3*one(t), zero(t)])
    @test t*xHt == HomogeneousPolynomial([t, zero(t)])
    @test complex(0,1)*xHt == HomogeneousPolynomial([1im*one(t), zero(1im*t)])
    @test eltype(complex(0,1)*xHt) == HomogeneousPolynomial{Taylor1{Complex{Float64}}}
    @test TS.numtype(complex(0,1)*xHt) == Taylor1{Complex{Float64}}
    @test (xHt+yHt)(1, 1) == 1+t
    @test (xHt+yHt)([1, 1]) == (xHt+yHt)((1, 1))

    tN1 = TaylorN([HomogeneousPolynomial([t]), xHt, yHt^2])
    @test findfirst(tN1) == 0
    @test tN1[0] == HomogeneousPolynomial([t])
    @test tN1(t,one(t)) == 2t+t^2
    @test findfirst(tN1(t,one(t))) == 1
    @test tN1([t,one(t)]) == tN1((t,one(t)))
    t1N = convert(Taylor1{TaylorN{Float64}}, tN1)
    @test findfirst(zero(tN1)) == -1
    @test t1N[0] == HomogeneousPolynomial(1)
    ctN1 = convert(TaylorN{Taylor1{Float64}}, t1N)
    @test convert(eltype(tN1), tN1) === tN1
    @test eltype(tN1) == TaylorN{Taylor1{Float64}}
    @test eltype(Taylor1([xH])) == Taylor1{HomogeneousPolynomial{Int}}
    @test TS.numtype(xHt) == Taylor1{Float64}
    @test TS.numtype(tN1) == Taylor1{Float64}
    @test TS.numtype(Taylor1([xH])) == HomogeneousPolynomial{Int}
    @test TS.numtype(t1N) == TaylorN{Float64}
    @test normalize_taylor(tN1) == tN1
    @test order(HomogeneousPolynomial([Taylor1(1), 1.0+Taylor1(2)])) == 1
    @test 3*tN1 == TaylorN([HomogeneousPolynomial([3t]),3xHt,3yHt^2])
    @test t*tN1 == TaylorN([HomogeneousPolynomial([t^2]),xHt*t,t*yHt^2])
    @test string(tN1) ==
        " ( 1.0 t + 𝒪(t⁴)) + ( 1.0 + 𝒪(t⁴)) x₁ + ( 1.0 t² + 𝒪(t⁴)) x₂² + 𝒪(‖x‖³)"
    @test string(t1N) ==
        "  1.0 x₁ + 𝒪(‖x‖³) + ( 1.0 + 𝒪(‖x‖³)) t + ( 1.0 x₂² + 𝒪(‖x‖³)) t² + 𝒪(t⁴)"
    @test tN1 == ctN1
    @test tN1+tN1 == 2*tN1
    @test tN1+1im*tN1 == complex(1,1)*tN1
    @test tN1+t == t+tN1
    @test tN1-t == -t+tN1
    zeroN1 = zero(tN1)
    oneN1 = one(tN1)
    @test tN1-tN1 == zeroN1
    @test zero(zeroN1) == zeroN1
    @test zeroN1 == zero.(tN1)
    @test oneN1[0] == one(tN1[0])
    @test oneN1[1:end] == zero.(tN1[1:end])
    for i in eachindex(tN1.coeffs)
        @test order(tN1.coeffs[i]) == order(zeroN1.coeffs[i]) == order(oneN1.coeffs[i])
    end
    @test string(t1N*t1N) ==
        "  1.0 x₁² + 𝒪(‖x‖³) + ( 2.0 x₁ + 𝒪(‖x‖³)) t + ( 1.0 + 𝒪(‖x‖³)) t² + ( 2.0 x₂² + 𝒪(‖x‖³)) t³ + 𝒪(t⁴)"
    @test !(@inferred isnan(tN1))
    @test !(@inferred isinf(tN1))

    @testset "Lexicographic order: HomogeneousPolynomial{Taylor1{T}} and TaylorN{Taylor1{T}}" begin
        @test 1 > xHt > yHt > xHt^2 > 0
        @test -xHt^2 < -xHt*yHt < -yHt^2 < -xHt^3 < -yHt^3 < 0
        @test 1 ≥ tN1 > 2*tN1^2 > 100*tN1^3 > 0
        @test -2*tN1 < -tN1^2 ≤ 0
    end

    @test round(tN1 + 1.0, RoundDown) === round(t1N + 1.0, RoundDown) === round(1.0, RoundDown)
    @test round(Int, tN1 + 1.0) === round(Int, t1N + 1.0) === round(Int, 1.0)
    @test ceil(tN1 + 1.0) === ceil(t1N + 1.0) === ceil(1.0)
    @test ceil(Int, tN1 + 1.0) === ceil(Int, t1N + 1.0) === ceil(Int, 1.0)
    @test floor(tN1 + 1.0) === floor(t1N + 1.0) === floor(1.0)
    @test floor(Int, tN1 + 1.0) === floor(Int, t1N + 1.0) === floor(Int, 1.0)
    @test mod(tN1+1,1.0) == 0+tN1
    @test mod(tN1-1.125,2) == 0.875+tN1
    @test (rem(tN1+1.125,1.0))[0][1] == 0.125 + t
    @test (rem(tN1-1.125,2))[0][1] == -1.125 + t
    @test mod2pi(-3pi+tN1)[0][1][0] ≈ pi
    @test mod2pi(0.125+2pi+tN1)[0][1][0] ≈ 0.125
    @test mod(t1N+1.125,1.0) == 0.125+t1N
    @test mod(t1N-1.125,2) == 0.875+t1N
    @test (rem(t1N+1.125,1.0))[0] == 0.125 + t1N[0]
    @test (rem(t1N-1.125,2))[0] == -1.125 + t1N[0]
    @test mod2pi(-3pi+t1N)[0][0][1] ≈ pi
    @test mod2pi(0.125+2pi+t1N)[0][0][1] ≈ 0.125

    @test abs(tN1+1) == 1+tN1
    @test abs(tN1-1) == 1-tN1
    @test_throws DomainError abs(tN1)
    @test_throws DomainError abs(t1N)

    @test abs2(im*(tN1+1)) == (1+tN1)^2
    @test abs2(im*(tN1-1)) == (1-tN1)^2
    @test abs(im*(tN1+1)) == 1+tN1
    @test abs(im*(tN1-1)) == 1-tN1

    @test convert(Array{Taylor1{TaylorN{Float64}},1}, [tN1, tN1]) == [t1N, t1N]
    @test convert(Array{Taylor1{TaylorN{Float64}},2}, [tN1 tN1]) == [t1N t1N]
    @test convert(Array{TaylorN{Taylor1{Float64}},1}, [t1N, t1N]) == [tN1, tN1]
    @test convert(Array{TaylorN{Taylor1{Float64}},2}, [t1N t1N]) == [tN1 tN1]

    @test evaluate(t1N, 0.0) == TaylorN(xH, 2)
    @test t1N() == TaylorN(xH, 2)
    @test string(evaluate(t1N, 0.0)) == " 1.0 x₁ + 𝒪(‖x‖³)"
    @test string(evaluate(t1N^2, 1.0)) == " 1.0 + 2.0 x₁ + 1.0 x₁² + 2.0 x₂² + 𝒪(‖x‖³)"
    @test string((t1N^(2//1))(1.0)) == " 1.0 + 2.0 x₁ + 1.0 x₁² + 2.0 x₂² + 𝒪(‖x‖³)"
    v = TaylorN.(zeros(Float64, 2), 2)
    vt1N = [t1N, t1N^2]
    @test eltype(evaluate(vt1N, 0.5im)) == TaylorN{ComplexF64}
    @test isnothing(evaluate!(vt1N, 0.0, v))
    @test v == [TaylorN(1), TaylorN(1)^2]
    evaluate!(vt1N, 0.5, v)
    @test v == evaluate(vt1N, 0.5)
    @test (@allocated evaluate!(vt1N, 0.5, v)) == 0
    δtN = 0.5 + t1N[0]
    auxN = zero(v[1])
    evaluate!(vt1N, δtN, v)
    @test v == evaluate(vt1N, δtN)
    evaluate!(vt1N, δtN, v, auxN)
    @test v == evaluate(vt1N, δtN)
    @test (@allocated evaluate!(vt1N, δtN, v, auxN)) == 0
    mixedN_constant = TaylorN(t1N[0].space, 1.0, 0)
    mixedN = Taylor1([mixedN_constant, t1N[0]], 1)
    evaluate!(mixedN, 1 + t1N[0], v[1], auxN)
    @test v[1] == evaluate(mixedN, 1 + t1N[0])
    limitedN_coeff = TaylorN(t1N[0].space, 1.0, 1)
    limitedN = Taylor1([limitedN_coeff, t1N[0]], 1)
    @test_throws DimensionMismatch evaluate!(limitedN, 0.5, v[1])
    @test_throws DimensionMismatch evaluate!(limitedN, 1 + t1N[0], v[1], auxN)
    @test tN1() == t
    @test evaluate(tN1, :x₁ => 1.0) == TaylorN([HomogeneousPolynomial([1.0+t]), zero(xHt), yHt^2])
    @test evaluate(tN1, 1, 1.0) == TaylorN([HomogeneousPolynomial([1.0+t]), zero(xHt), yHt^2])
    @test evaluate(t, t1N) == t1N
    @test evaluate(t1N, 0.5) == t1N[0] + t1N[1]/2 + t1N[2]/4
    @test evaluate(t1N, [t1N[0], zero(t1N[0])]) == Taylor1([t1N[0], t1N[1]], order(t))
    @test evaluate(t1N, 2, 0.0) == Taylor1([t1N[0], t1N[1]], order(t))
    @test evaluate(t1N, 1, 0.0) == Taylor1([zero(t1N[0]), t1N[1], t1N[2]], order(t))

    # Tests for functions of mixtures
    t1N = Taylor1([zero(TaylorN(Float64,1)), one(TaylorN(Float64,1))], 6)
    t = Taylor1(3)
    xHt = HomogeneousPolynomial([one(t), zero(t)])
    yHt = HomogeneousPolynomial([zero(t), t])
    x = TaylorN(1, order=2)
    y = TaylorN(2, order=2)
    xN1 = TaylorN([HomogeneousPolynomial(zero(t), 0), xHt, zero(yHt)], 2)
    yN1 = TaylorN([HomogeneousPolynomial(zero(t), 0), zero(xHt), yHt], 2)
    for fn in (exp, log, log1p, sin, cos, tan, sinh, cosh, tanh, asin, acos, atan, asinh, acosh, atanh)
        if fn == asin || fn == acos || fn == atanh || fn == log1p
            cc = 0.5
        elseif fn == asinh || fn == acosh
            cc = 1.5
        else
            cc = 1.0
        end
        @test x*fn(cc+t1N) == fn(cc+t)*xN1
        @test t*fn(cc+xN1) == fn(cc+x)*t1N
    end
    ee = zero(Taylor1(t1N[0:5], 6))
    ee_orig = exp(t1N)
    for ord in eachindex(t1N)
        TS.differentiate!(ee, ee_orig, ord)
    end
    @test iszero(ee[6])
    @test getcoeff.(ee, 0:5) == getcoeff.(ee_orig, 0:5)
    ee = differentiate(t1N, order(t1N))
    @test iszero(ee)
    @test iszero(order(ee))

    vt = zeros(Taylor1{Float64},2)
    @test isnothing(evaluate!([tN1, tN1^2], [t, t], vt))
    @test vt == [2t, 4t^2]

    tint = Taylor1(Int, 10)
    t = Taylor1(10)
    x = TaylorN( [HomogeneousPolynomial(zero(t), 5), HomogeneousPolynomial([one(t),zero(t)])], 5)
    y = TaylorN(typeof(tint), 2, order=5)
    @test typeof(x) == TaylorN{Taylor1{Float64}}
    @test eltype(y) == TaylorN{Taylor1{Int}}
    @test TS.numtype(y) == Taylor1{Int}
    @test -x == 0 - x
    @test +y == y
    @test one(y)/(1+x) == 1 - x + x^2 - x^3 + x^4 - x^5
    @test one(y)/(1+y) == 1 - y + y^2 - y^3 + y^4 - y^5
    @test (1+y)/one(t) == 1 + y
    @test typeof(y+t) == TaylorN{Taylor1{Float64}}

    t = Taylor1(4)
    xN, yN = variables()
    @test evaluate(1.0 + t + t^2, xN) == 1.0 + xN + xN^2
    v1 = [1.0 + t + t^2 + t^4, 1.0 - t^2 + t^3]
    @test v1(yN^2) == [1.0 + yN^2 + yN^4, 1.0 - yN^4 + yN^6]
    tN = Taylor1([zero(xN), one(xN)], 4)
    q1N = 1 + yN*tN + xN*tN^4
    @test q1N(-1.0) == 1.0 - yN + xN
    @test q1N(-xN^2) == 1.0 - xN^2*yN

    # See #92 and #94
    δx, δy = variables!("δx δy", nowarn=true)
    xx = 1+Taylor1(δx, 5)
    yy = 1+Taylor1(δy, 5)
    tt = Taylor1([zero(δx), one(δx)], order(xx))
    @test all((xx, yy) .== (TS.fixorder(xx, yy)))
    @test typeof(xx) == Taylor1{TaylorN{Float64}}
    @test eltype(xx) == Taylor1{TaylorN{Float64}}
    @test TS.numtype(tt) == TaylorN{Float64}
    @test !(@inferred isnan(xx))
    @test !(@inferred isnan(δx))
    @test !(@inferred isinf(xx))
    @test !(@inferred isinf(δx))
    @test +xx == xx
    @test -xx == 0 - xx
    @test xx/1.0 == 1.0*xx
    @test xx + xx == xx*2
    @test xx - xx == zero(xx)
    @test xx*xx == xx^2
    @test xx/xx == one(xx)
    @test xx*δx + Taylor1(typeof(δx),5) == δx + δx^2 + Taylor1(typeof(δx),5)
    @test xx/(1+δx) == one(xx)
    @test 1/(1-tt) == 1 + tt + tt^2 + tt^3 + tt^4 + tt^5
    @test xx/(1-tt) == xx * (1 + tt + tt^2 + tt^3 + tt^4 + tt^5)
    res = xx * tt
    @test 1/(1-xx*tt) == 1 + res + res^2 + res^3 + res^4 + res^5
    @test typeof(xx+δx) == Taylor1{TaylorN{Float64}}
    res = 1/(1+δx)
    @test one(xx)/(xx*(1+tt)) == Taylor1([res, -res, res, -res, res, -res])
    res = 1/(1+δx)^2
    @test (xx^2 + yy^2)/(xx*yy) == xx/yy + yy/xx
    @test ((xx+yy)*tt)^2/((xx+yy)*tt) == (xx+yy)*tt
    @test sqrt(xx) == Taylor1(sqrt(xx[0]), order(xx))
    @test xx^0.25 == sqrt(sqrt(xx))
    @test (xx*yy*tt^2)^0.5 == sqrt(xx*yy)*tt
    FF(x,y,t) = (1 + x + y + t)^4
    QQ(x,y,t) = x^4.0 + (4*y + (4*t + 4))*x^3 + (6*y^2 + (12*t + 12)*y + (6*t^2 + 12*t + 6))*x^2 +
        (4*y^3 + (12*t + 12)*y^2 + (12*t^2 + 24*t + 12)*y + (4*t^3 + 12*t^2 + 12*t + 4))*x +
        (y^4 + (4*t + 4)*y^3 + (6*t^2 + 12*t + 6)*y^2 + (4*t^3 + 12*t^2 + 12*t + 4)*y +
            (t^4 + 4*t^3 + 6*t^2 + 4*t + 1))
    @test FF(xx, yy, tt) == QQ(xx, yy, tt)
    @test FF(tt, yy-1, xx-1) == QQ(xx-1, yy-1, tt)
    @test (xx+tt)^4 == xx^4 + 4*xx^3*tt + 6*xx^2*tt^2 + 4*xx*tt^3 + tt^4
    pp = xx*yy*(1+tt)^4
    @test pp^0.25 == sqrt(sqrt(pp))
    @test (xx*yy)^(3/2)*(1+tt+tt^2) == (sqrt(xx*yy))^3*(1+tt+tt^2)
    @test sqrt((xx+yy+tt)^3) ≈ (xx+yy+tt)^1.5
    @test (sqrt(xx*(1+tt)))^4 ≈ xx^2 * (1+tt)^2
    @test (sqrt(xx*(1+tt)))^5 ≈ xx^2.5 * (1+tt)^2.5
    @test exp(xx) == exp(xx[0]) + zero(tt)
    @test exp(xx+tt) ≈ exp(xx)*exp(tt)
    @test norm(exp(xx+tt) - exp(xx)*exp(tt), Inf) < 5e-17
    @test expm1(xx) == expm1(xx[0]) + zero(tt)
    @test expm1(xx+tt) ≈ exp(xx+tt)-1
    @test norm(expm1(xx+tt) - (exp(xx+tt)-1), Inf) < 5e-16
    @test log(xx) == log(xx[0]) + zero(tt)
    @test log((xx+tt)*(yy+tt)) ≈ log(xx+tt)+log(yy+tt)
    @test log(pp) ≈ log(xx)+log(yy)+4*log(1+tt)
    @test norm(log(pp) - (log(xx)+log(yy)+4*log(1+tt)), Inf) < 1e-15
    @test log1p(xx) == log1p(xx[0]) + zero(tt)
    @test log1p(xx+tt) == log(1+xx+tt)
    exp(log(xx+tt)) == log(exp(xx+tt))
    @test exp(log(pp)) ≈ log(exp(pp))
    @test sincos(δx + xx) == sincos(δx+xx[0]) .+ zero(tt)
    qq = sincos(pp)
    @test exp(im*pp) ≈ qq[2] + im*qq[1]
    @test qq[2]^2 ≈ 1 - qq[1]^2
    @test sincospi(δx + xx - 1) == sincospi(δx + xx[0] - 1) .+ zero(tt)
    @test all(sincospi(pp) .≈ sincos(pi*pp))
    @test tan(xx) == tan(xx[0]) + zero(tt)
    @test tan(xx+tt) ≈ sin(xx+tt)/cos(xx+tt)
    @test tan(xx*tt) ≈ sin(xx*tt)/cos(xx*tt)
    @test asin(xx-1) == asin(xx[0]-1) + zero(tt)
    @test asin(sin(xx+tt)) ≈ xx + tt
    @test sin(asin(xx*tt)) == xx * tt
    @test_throws DomainError asin(2+xx+tt)
    @test acos(xx-1) == acos(xx[0]-1) + zero(tt)
    @test acos(cos(xx+tt)) ≈ xx + tt
    @test cos(acos(xx*tt)) ≈ xx * tt
    @test_throws DomainError acos(2+xx+tt)
    @test atan(xx) == atan(xx[0]) + zero(tt)
    @test atan(tan(xx+tt)) ≈ xx + tt
    @test tan(atan(xx*tt)) == xx * tt
    @test TS.sinhcosh(xx) == TS.sinhcosh(xx[0]) .+ zero(tt)
    qq = TS.sinhcosh(pp)
    @test qq[2]^2 - 1 ≈ qq[1]^2
    @test 2*qq[1] ≈ exp(pp)-exp(-pp)
    @test 2*qq[2] ≈ exp(pp)+exp(-pp)
    qq = TS.sinhcosh(xx+tt)
    @test tanh(xx) == tanh(xx[0]) + zero(tt)
    @test tanh(xx+tt) ≈ qq[1]/qq[2]
    qq = TS.sinhcosh(xx*tt)
    @test tanh(xx*tt) ≈ qq[1]/qq[2]
    @test asinh(xx-1) == asinh(xx[0]-1) + zero(tt)
    @test asinh(sinh(xx+tt)) ≈ xx + tt
    @test sinh(asinh(xx*tt)) == xx * tt
    @test acosh(xx+1) == acosh(xx[0]+1) + zero(tt)
    @test acosh(cosh(xx+1+tt)) ≈ xx + 1 + tt
    @test cosh(acosh(2+xx*tt)) ≈ 2 + xx * tt
    @test_throws DomainError acosh(xx+tt)
    @test atanh(xx-1) == atanh(xx[0]-1) + zero(tt)
    @test atanh(tanh(-1+xx+tt)) ≈ -1 + xx + tt
    @test tanh(atanh(xx*tt)) ≈ xx * tt
    @test_throws DomainError atanh(xx+tt)

    # pp = xx*yy*(1+tt)^4
    @test evaluate(pp, 1, 0.0) == yy*(1+tt)^4
    @test evaluate(pp, 2, 0.0) == xx*(1+tt)^4
    @test evaluate(t, tt) == tt
    @test evaluate(tt, t) == tt
    @test evaluate(xx, 2, δy) == xx
    @test evaluate(xx, 1, δy) == yy

    #testing evaluate and function-like behavior of Taylor1, TaylorN for mixtures:
    t = Taylor1(25)
    p = cos(t)
    q = sin(t)
    a = [p,q]
    dx = variables!("x", numvars=4, order=10, nowarn=true)
    P = sin.(dx)
    v = [1.0,2,3,4]
    F(x) = [sin(sin(x[4]+x[3])), sin(cos(x[3]-x[2])), cos(sin(x[1]^2+x[2]^2)), cos(cos(x[2]*x[3]))]
    Q = F(v+dx)
    diff_evals = cos(sin(dx[1]))-p(P[1])
    @test norm(diff_evals, Inf) < 1e-15
    #evaluate a Taylor1 at a TaylorN
    @test p(P) == evaluate(p, P)
    @test q(Q) == evaluate(q, Q)
    #evaluate an array of Taylor1s at a TaylorN
    aT1 = [p,q,p^2,log(1+q)] #an array of Taylor1s
    @test aT1(Q[4]) == evaluate(aT1, Q[4])
    @test (aT1.^2)(Q[3]) == evaluate(aT1.^2, Q[3])
    #evaluate a TaylorN at an array of Taylor1s
    @test P[1](aT1) == evaluate(P[1], aT1)
    @test P[1](aT1) == evaluate(P[1], (aT1...,))
    @test Q[2](aT1) == evaluate(Q[2], [aT1...])
    #evaluate an array of TaylorN{Float64} at an array of Taylor1{Float64}
    @test P(aT1) == evaluate(P, aT1)
    @test Q(aT1) == evaluate(Q, aT1)
    #test evaluation of an Array{TaylorN{Taylor1}} at an Array{Taylor1}
    aH1 = [
        HomogeneousPolynomial([Taylor1(rand(2))]),
        HomogeneousPolynomial([Taylor1(rand(2)),Taylor1(rand(2)),
            Taylor1(rand(2)),Taylor1(rand(2))])
        ]
    bH1 = [
        HomogeneousPolynomial([Taylor1(rand(2))]),
        HomogeneousPolynomial([Taylor1(rand(2)),Taylor1(rand(2)),
            Taylor1(rand(2)),Taylor1(rand(2))])
        ]
    aTN1 = TaylorN(aH1); bTN1 = TaylorN(bH1)
    x = [aTN1, bTN1]
    δx = [Taylor1(rand(3)) for i in 1:4]
    @test typeof(x) == Array{TaylorN{Taylor1{Float64}},1}
    @test typeof(δx) == Array{Taylor1{Float64},1}
    x0 = Array{Taylor1{Float64}}(undef, length(x))
    eval_x_δx = evaluate(x,δx)
    @test x(δx) == eval_x_δx
    evaluate!(x,δx,x0)
    @test x0 == eval_x_δx
    δx_num = rand(4)
    eval_x_δx_num = evaluate.(x, Ref(δx_num))
    @test evaluate(x, δx_num) == eval_x_δx_num
    @test evaluate(x, (δx_num...,)) == eval_x_δx_num
    x0_num = similar(eval_x_δx_num)
    evaluate!(x, δx_num, x0_num)
    @test x0_num == eval_x_δx_num
    @test typeof(evaluate(x[1],δx)) == Taylor1{Float64}
    @test x() == map(y->y[0][1], x)
    for i in eachindex(x)
        @test evaluate(x[i],δx) == eval_x_δx[i]
        @test x[i](δx) == eval_x_δx[i]
    end
    p11 = Taylor1([sin(t),cos(t)])
    @test evaluate(p11,t) == sin(t)+t*cos(t)
    @test p11(t) == sin(t)+t*cos(t)
    a11 = Taylor1([t,t^2,exp(-t),sin(t),cos(t)])
    b11 = t+t*(t^2)+(t^2)*(exp(-t))+(t^3)*sin(t)+(t^4)*cos(t)
    diff_a11b11 = a11(t)-b11
    @test norm(diff_a11b11.coeffs, Inf) < 1E-19

    X, Y = variables!(Taylor1{Float64}, "x y", nowarn=true)
    @test typeof( norm(X) ) == Float64
    @test norm(X) > 0
    @test norm(X+Y) == sqrt(2)
    @test norm(-10X+4Y,Inf) == 10.


    X,Y = convert(Taylor1{TaylorN{Float64}},X), convert(Taylor1{TaylorN{Float64}},Y)
    @test typeof( norm(X) ) == Float64
    @test norm(X) > 0
    @test norm(X+Y) == sqrt(2)
    @test norm(-10X+4Y,Inf) == 10.


    @test TS.rtoldefault(TaylorN{Taylor1{Int}}) == 0
    @test TS.rtoldefault(Taylor1{TaylorN{Int}}) == 0
    for T in (Float64, BigFloat)
        @test TS.rtoldefault(TaylorN{Taylor1{T}}) == sqrt(eps(T))
        @test TS.rtoldefault(Taylor1{TaylorN{T}}) == sqrt(eps(T))
        @test TS.real(TaylorN{Taylor1{T}}) == TaylorN{Taylor1{T}}
        @test TS.real(Taylor1{TaylorN{T}}) == Taylor1{TaylorN{T}}
        @test TS.real(TaylorN{Taylor1{Complex{T}}}) == TaylorN{Taylor1{T}}
        @test TS.real(Taylor1{TaylorN{Complex{T}}}) == Taylor1{TaylorN{T}}
    end

    rndT1(ord1) = Taylor1(-1 .+ 2rand(ord1+1)) # generates a random Taylor1 with order `ord`
    nmonod(s, d) = binomial(d+s-1, d) #number of monomials in s variables with exact degree d
    #rndHP generates a random `ordHP`-th order homog. pol. of Taylor1s, each with order `ord1`
    rndHP(ordHP, ord1) = HomogeneousPolynomial( [rndT1(ord1) for i in 1:nmonod(get_numvars(), ordHP)] )
    #rndTN generates a random `ordHP`-th order TaylorN of of Taylor1s, each with order `ord1`
    rndTN(ordN, ord1) = TaylorN([rndHP(i, ord1) for i in 0:ordN])

    P = rndTN(order(), 3)
    @test P ≈ P
    Q = deepcopy(P)
    Q[2][2] = Taylor1([NaN, Inf])
    @test (@inferred isnan(Q))
    @test (@inferred isinf(Q))
    @test !isfinite(Q)
    Q[2][2] = P[2][2]+sqrt(eps())/2
    @test isapprox(P, Q, rtol=1.0)
    Q[2][2] = P[2][2]+10sqrt(eps())
    @test !isapprox(P, Q, atol=sqrt(eps()), rtol=0)
    @test P ≉ Q^2
    Q[2][2] = P[2][2]+eps()/2
    @test isapprox(Q, Q, atol=eps(), rtol=0)
    @test isapprox(Q, P, atol=eps(), rtol=0)
    Q[2][1] = P[2][1]-10eps()
    @test !isapprox(Q, P, atol=eps(), rtol=0)
    @test P ≉ Q^2

    X, Y = variables!(BigFloat, "x y", numvars=2, order=6, nowarn=true)
    p1N = Taylor1([X^2,X*Y,Y+X,Y^2])
    q1N = Taylor1([X^2,(1.0+sqrt(eps(BigFloat)))*X*Y,Y+X,Y^2])
    @test p1N ≈ p1N
    @test p1N ≈ q1N

    Pv = [rndTN(order(), 3), rndTN(order(), 3)]
    Qv = convert.(Taylor1{TaylorN{Float64}}, Pv)

    @test TS.jacobian(Pv) == TS.jacobian(Qv)

    @test_throws ArgumentError Taylor1(2) + TaylorN(1)
    @test_throws ArgumentError Taylor1(2) - TaylorN(1)
    @test_throws ArgumentError Taylor1(2) * TaylorN(1)
    @test_throws ArgumentError TaylorN(2) / Taylor1(1)

    # Issue #342 and PR #343
    z0N = -1.333+variables()[1]
    z = Taylor1(z0N,20)
    z[20][1][1] = 5.0
    @test z[0][0][1] == -1.333
    @test z[20][1][1] == 5.0
    for i in 1:19
        for j in eachindex(z[i].coeffs)
            @test all(z[i].coeffs[j][:] .== 0.0)
        end
    end
    @test all(z[20][1][2:end] .== 0.0)
    intz = integrate(z)
    intz[20] = z[0]
    @test intz[1] == z[0]
    @test intz[20] == z[0]
    for i in 2:19
        @test iszero(intz[i])
    end

    a = sum(exp.(variables()).^2)
    b = Taylor1([a])
    bcopy = deepcopy(b)
    c = Taylor1(constant_term(b),0)
    c[0][0][1] = 0.0
    b == bcopy

    #347
    a = Taylor1([1.0+X,-X, Y, X-Y,X])
    b = deepcopy(a)
    b[0] = zero(a[0])
    b.coeffs[2:end] .= zero(b.coeffs[1])
    @test iszero(b)
    @test b.coeffs[2] === b.coeffs[3]
    b.coeffs[2:end] .= zero.(b.coeffs[1])
    @test !(b.coeffs[2] === b.coeffs[3])
    x = Taylor1([1.0+X,-X, Y, X-Y,X])
    z = zero(x)
    two = 2one(x[0])
    @test two/x == 2/x == 2.0/x
    @test (2one(x))/x == 2/x

    dq = variables()
    x = Taylor1(exp.(dq), 5)
    x[1] = sin(dq[1]*dq[2])
    @test x[1] == sin(dq[1]*dq[2])
    @test x[1] !== sin(dq[1]*dq[2])

    @testset "Promotion of Taylor1 with Taylor1{TaylorN} in a non-default JetSpace" begin
        sploc = JetSpace(4, ["x"])   # ADAPT: any space different from TS.default_space
        ordN = 4
        x  = Taylor1([TaylorN(sploc, float(k+1), ordN) for k in 0:4]) # order 5
        t  = Taylor1([1.0, 2.0, 3.0, 4.0])                            # order 3
        ti = Taylor1([1, 2, 3, 4])                                    # Int coeffs
        tl = Taylor1([TaylorN(sploc, c, ordN) for c in t.coeffs])     # t lifted by hand
        @test sploc != TS.default_space[]

        # promote itself
        a, b = promote(t, x)
        @test a isa Taylor1{TaylorN{Float64}} && b isa Taylor1{TaylorN{Float64}}
        @test all(space(a[k]) == sploc for k in eachindex(a))
        @test reverse(promote(x, t)) == (a, b)

        # the original bug
        tN = t + zero(x)
        @test all(space(tN[k]) == sploc for k in eachindex(tN))
        @test tN == tl + zero(x)
        for op in (+, -), u in (t, ti)
            @test op(u, x) == op(tl, x)
            @test op(x, u) == op(x, tl)
            @test order(op(u, x)) == order(u)
        end
        @test t * x == tl * x
        @test x / (1 + t) == x / (1 + tl)
    end

    @testset "Test Base.float overloads for Taylor1 and TaylorN mixtures" begin
        q = variables(Int)
        x1N = Taylor1(q)
        @test float(x1N) == Taylor1(float.(q))
        xN1 = convert(TaylorN{Taylor1{Int}}, x1N)
        @test float(xN1) == convert(TaylorN{Taylor1{Float64}}, Taylor1(float.(q)))
        @test float(Taylor1{TaylorN{Int}}) == Taylor1{TaylorN{Float64}}
        @test float(TaylorN{Taylor1{Int}}) == TaylorN{Taylor1{Float64}}
        @test float(TaylorN{Taylor1{Complex{Int}}}) == TaylorN{Taylor1{Complex{Float64}}}
    end
end

@testset "Tests with nested Taylor1s" begin
    set_taylor1_varname(2, " t s")
    @test TS._params_Taylor1_.num_vars == 2
    @test TS._params_Taylor1_.var_name == ["t", "s"]
    set_taylor1_varname(2, [" t", "s w"])
    @test TS._params_Taylor1_.num_vars == 2
    @test TS._params_Taylor1_.var_name == ["t", "s"]
    ti = Taylor1(3)
    tii = Taylor1([zero(ti), one(ti)], 9)
    @test typeof(promote(ti, tii)[1]) == Taylor1{Taylor1{Float64}}
    @test findfirst(tii) == 1
    @test TS.numtype(tii) == Taylor1{Float64}
    @test normalize_taylor(tii) == tii
    @test normalize_taylor(Taylor1([zero(tii), one(tii)], 5)) == Taylor1([zero(tii), one(tii)], 5)
    @test convert(eltype(tii), tii) === tii
    @test string(ti) == " 1.0 t + 𝒪(t⁴)"
    @test string(tii) == " ( 1.0 + 𝒪(t⁴)) s + 𝒪(s¹⁰)"
    @test string(tii^2) == " ( 1.0 + 𝒪(t⁴)) s² + 𝒪(s¹⁰)"
    @test ti + tii == Taylor1([ti, one(ti)], 9)
    titii = ti * tii
    @test string(titii) == " ( 1.0 t + 𝒪(t⁴)) s + 𝒪(s¹⁰)"
    @test titii == Taylor1([zero(ti), ti], 9)
    vtii = [tii, titii]
    @test eltype(evaluate(vtii, 0.5im)) == Taylor1{ComplexF64}
    vtii_dest = [zero(ti), zero(ti)]
    evaluate!(vtii, 0.5, vtii_dest)
    @test vtii_dest == evaluate(vtii, 0.5)
    @test (@allocated evaluate!(vtii, 0.5, vtii_dest)) == 0
    @testset "Numeric evaluation with deeper Taylor1 nesting" begin
        # Three and four Taylor1 layers use the general allocating fallback.
        for inner in (tii, Taylor1([tii, one(tii)], 1))
            src = [Taylor1([inner, 2inner], 1),
                   Taylor1([one(inner), inner], 1)]
            original = deepcopy(src)
            dest = [zero(inner), zero(inner)]
            @test isnothing(evaluate!(src, 0.5, dest))
            @test dest == [2inner, one(inner) + 0.5inner]
            @test isnothing(evaluate!(src, 2.0, view(dest, :)))
            @test dest == [5inner, one(inner) + 2inner]
            @test src == original
        end
    end
    δt1 = 0.5 + ti
    aux1 = zero(vtii_dest[1])
    evaluate!(vtii, δt1, vtii_dest)
    @test vtii_dest == evaluate(vtii, δt1)
    evaluate!(vtii, δt1, vtii_dest, aux1)
    @test vtii_dest == evaluate(vtii, δt1)
    @test (@allocated evaluate!(vtii, δt1, vtii_dest, aux1)) == 0
    @test_throws DimensionMismatch evaluate!(vtii, δt1, Taylor1{Float64}[])
    mixed_inner_orders = Taylor1([Taylor1(1.0, 0), ti], 1)
    evaluate!(mixed_inner_orders, 1 + ti, vtii_dest[1], aux1)
    @test vtii_dest[1] == evaluate(mixed_inner_orders, 1 + ti)
    constant_inner_coeffs = Taylor1([Taylor1(1.0, 0), Taylor1(2.0, 0)], 1)
    inner_variable = Taylor1(4)
    @test evaluate(constant_inner_coeffs, 1 + inner_variable) ==
        3 + 2inner_variable
    limited_inner_order = Taylor1([Taylor1([1.0, 1.0], 1), ti], 1)
    @test_throws DimensionMismatch evaluate!(limited_inner_order, 1 + ti,
        vtii_dest[1], aux1)
    aliased_aux = 1 + ti
    aliased_source = Taylor1([one(ti), zero(ti)], 1)
    aliased_source.coeffs[2] = aliased_aux   # direct write: the constructor no longer keeps references
    @test_throws ArgumentError evaluate!(aliased_source, 1 + ti,
        zero(ti), aliased_aux)
    mixed_dest = Taylor1(0.0, 0)
    mixed_aux = zero(mixed_dest)
    evaluate!(mixed_inner_orders, 1 + ti, mixed_dest, mixed_aux)
    @test mixed_dest[0] == evaluate(mixed_inner_orders, 1 + ti)[0]
    @test titii / tii == ti
    @test order(titii/tii) == order(tii)-1
    @test titii / ti == tii
    @test order(titii/ti) == order(tii)
    @test ti^2-tii^2 == (ti+tii)*(ti-tii)
    @test findfirst(ti^2-tii^2) == 0
    @test sin(tii) ≈ Taylor1(one(ti) .* sin(Taylor1(10)).coeffs, 9)
    @test tii(1 + ti) == 1 + ti
    @test tii(1 + ti) isa Taylor1{Float64}
    @test ti(1 + tii) == 1 + tii
    @test constant_term(ti+tii) == ti
    @test linear_polynomial(ti*tii) == Taylor1([zero(ti), ti], 9)
    @test order(linear_polynomial(tii)) == order(tii)
    @test nonlinear_polynomial(tii+ti*tii^2) == Taylor1([zero(ti), zero(ti), ti], 9)
    @test ti(1 + tii) isa Taylor1{Taylor1{Float64}}
    @test sqrt(titii^2) == titii
    @test order(sqrt(titii^2)) == order(tii) >> 1
    @test (titii^3)^(1/3) == titii
    @test order(sqrt(titii^2)) == order(tii) >> 1
    ti2tii = ti^2 * tii
    tti = (ti2tii/tii)/ti
    @test order(tti) == order(tii)-1
    @test order(tti[0]) == order(ti)-1
    @test isapprox(abs2(exp(im*tii)), one(tii))
    @test isapprox(abs(exp(im*tii)), one(tii))
    tii = Taylor1([1/(1+ti), one(ti)], 9)
    @test tii(1.0) == 1 + 1/(1+ti)
    @test cos(tii)(0.0) == cos(tii[0])
    @test tii(ti) == tii[0] + ti
    @test evaluate(tii*ti, ti) == tii[0]*ti + ti^2

    # Testing automatic setting of TS._params_taylor1_
    tii = Taylor1([zero(ti), one(ti)], 9)
    tiii = Taylor1([zero(tii), one(tii)], 9)
    @test string(tiii) == " (  1.0 + 𝒪(t₁⁴) + 𝒪(t₂¹⁰)) t₃ + 𝒪(t₃¹⁰)"
    @test TS._params_Taylor1_.var_name == ["t₁", "t₂", "t₃"]
    @test TS._params_Taylor1_.num_vars == 3
    # The next tests are related to iss #326
    @test ti > ti^2 > tii > 0.0
    @test tiii < tii^2 < titii < ti^2

    @testset "Test setindex! method for nested Taylor1s" begin
        t = Taylor1(2)
        y = one(t)
        x = Taylor1([t,2t,t^2,0t,t^3])
        x[3] = y
        @test x[3] !== y
        y[2] = -5.0
        @test x[3][2] == 0.0
    end

    @testset "Tests for functions for nested Taylor1s" begin
        funs = (:exp, :expm1, :log, :log1p, :sin, :cos, :sinpi, :cospi, :tan, :asin, :acos, :atan,
            :sinh, :cosh, :tanh, :asinh, :atanh, :acosh)

        for fn in funs
            @eval begin
                c = 0.0
                ti = Taylor1(9)
                tii = Taylor1([zero(ti), one(ti)], 9)
                tiii = Taylor1([zero(tii), one(tii)], 9)
                t1 = ti + ti^2
                t2 = tii + tii^2
                t3 = tiii + tiii^2
                if $fn == log
                    c = 1.0
                elseif $fn == acosh
                    c = 2.0
                end
                for i in eachindex(t1)
                    ev_t1 = getcoeff($(fn)(c + t1), i)
                    ev_t2 = evaluate(getcoeff($(fn)(c + t2), i), 1.0)
                    ev_t3 = evaluate(getcoeff($(fn)(c + t3), i), 1.0)
                    @test ev_t1 ≈ ev_t2 == ev_t3
                end
                if $fn == log || $fn == acosh
                    c = 0.0
                end
            end
            # println()
        end

    end

    @testset "Tests for mutating functions of nested Taylor1s" begin
        import Random
        Random.seed!(1234)

        for (f, fc) in ((:+, :(add!)), (:-, :(subst!)), (:*, :(mul!)), (:/, :(div!)))
            @eval begin

                local inorder = 6
                local outorder = 25

                c = rand()
                x = Taylor1([Taylor1(rand(inorder+1), inorder) for _ in 1:outorder+1], outorder)
                y = Taylor1([Taylor1(rand(inorder+1), inorder) for _ in 1:outorder+1], outorder)
                z = zero(x)

                w = $f(x, y)
                for k in eachindex(z)
                    TS.$fc(z, x, y, k)
                end
                @test norm(z - w, Inf) == 0.0

                w = $f(c, y)
                for k in eachindex(z)
                    TS.$fc(z, c, y, k)
                end
                @test norm(z - w, Inf) == 0.0

                w = $f(y, c)
                for k in eachindex(z)
                    TS.$fc(z, y, c, k)
                end
                @test norm(z - w, Inf) == 0.0

            end
        end

        for r in (0, 1, 2, 0.5)
            local inorder = 6
            local outorder = 25

            x = Taylor1([Taylor1(rand(inorder+1), inorder) for _ in 1:outorder+1], outorder)
            y = zero(x)
            z = zero(x)

            w = x^r
            for k in eachindex(z)
                TS.pow!(z, x, y, r, k)
            end
            @test norm(z - w, Inf) == 0.0
        end
    end

    # Back to default
    set_taylor1_varname(1, "t")
    @test TS._params_Taylor1_.var_name == ["t"]
end

@testset "Test mixtures involving cis and cispi" begin
    sp = JetSpace(order = 4, variables = [:u, :v])
    u, v = variables(sp)
    t = Taylor1(8)
    real_series = (
        0.3 + t - 0.2 * t^2,                                            # Taylor1{Float64}
        0.3 + u - 0.5 * u * v + v^3,                                    # TaylorN{Float64}
        Taylor1([0.3 + u, 1.0 + 0.0 * u, 0.5 * v], 6),                  # Taylor1{TaylorN}
        Taylor1([Taylor1([0.3, 1.0], 3), Taylor1([0.0, 0.5], 3)], 6))   # Taylor1{Taylor1}
    for x in real_series
        z = cis(x)
        s, c = sincos(x)
        @test typeof(z) == typeof(im * x)
        @test real(z) == c
        @test imag(z) == s
        @test z ≈ exp(im * x)
        zpi = cispi(x)
        s, c = sincospi(x)
        @test typeof(zpi) == typeof(im * x)
        @test real(zpi) == c
        @test imag(zpi) == s
        @test zpi ≈ exp(im * pi * x)
    end

    complex_series = (
        (0.3 + 0.1im) + t - 0.2 * t^2,                                  # Taylor1{ComplexF64}
        (0.3 + 0.1im) + u - 0.5im * u * v)                              # TaylorN{ComplexF64}
    for x in complex_series
        s, c = sincos(x)
        @test cis(x) == c + im * s
        @test cis(x) ≈ exp(im * x)
        s, c = sincospi(x)
        @test cispi(x) == c + im * s
        @test cispi(x) ≈ exp(im * pi * x)
    end
end

# TaylorN iand Taylor1{Taylor1} represent the same expansions: the coefficient of
# TaylorN (ξ₁^i₁ ⋯ ξₙ^iₙ) is that of the nested Taylor1 whose innermost variable is
# the first one and outermost the last one, p[iₙ]⋯[i₁].
# The helpers below compare all coefficients that both expansions contain,
# up to rounding: relative tolerance `rtol`, and no absolute tolerance
# (coefficients that vanish in one representation must vanish in the other).
_coeff_like(p::TaylorN, idx::NTuple) = getcoeff(p, idx)
_coeff_like(p::Taylor1, idx::NTuple) = foldl((q, n) -> q[n], reverse(idx); init = p)
_orders_nested(p::Taylor1) = p[0] isa Taylor1 ? (_orders_nested(p[0])..., TS.order(p)) :
    (TS.order(p),)
function _same_coeffs(pN::TaylorN, pT::Taylor1; rtol = 32 * eps())
    ords = _orders_nested(pT)
    return all(isapprox(_coeff_like(pN, idx), _coeff_like(pT, idx); rtol)
        for idx in Iterators.product((0:n for n in ords)...) if sum(idx) ≤ TS.order(pN))
end
# Taylor1{TaylorN}: the first variables are those of the TaylorN coefficients, and the last one
# is that of the Taylor1; p[iₙ] is the TaylorN that holds the coefficients of ξₙ^iₙ
_coeff_like(p::Taylor1{<:TaylorN}, idx::NTuple) = getcoeff(p[idx[end]], idx[1:end-1])
function _same_coeffs(pN::TaylorN, pT::Taylor1{<:TaylorN}; rtol = 32 * eps())
    nv, M, K = get_numvars(pT[0]), TS.order(pT[0]), TS.order(pT)
    return all(isapprox(_coeff_like(pN, idx), _coeff_like(pT, idx); rtol)
        for idx in Iterators.product(ntuple(_ -> 0:M, nv)..., 0:K)
        if sum(idx[1:nv]) ≤ M && sum(idx) ≤ TS.order(pN))
end

@testset "Compare TaylorN and nested Taylor1s" begin
    # functions, and three expansion points within each domain
    _cmp_funs() = ((exp, (0.0, 0.3, -0.4)), (expm1, (0.0, 0.3, -0.4)), (log, (1.0, 1.5, 0.7)),
        (log1p, (0.0, 0.3, -0.4)), (sqrt, (1.0, 1.5, 0.7)), (sin, (0.0, 0.3, -0.4)),
        (cos, (0.0, 0.3, -0.4)), (tan, (0.0, 0.3, -0.4)), (sinpi, (0.0, 0.3, -0.4)),
        (cospi, (0.0, 0.3, -0.4)), (asin, (0.0, 0.3, -0.4)), (acos, (0.0, 0.3, -0.4)),
        (atan, (0.0, 0.3, -0.4)), (sinh, (0.0, 0.3, -0.4)), (cosh, (0.0, 0.3, -0.4)),
        (tanh, (0.0, 0.3, -0.4)), (asinh, (0.0, 0.3, -0.4)), (acosh, (2.0, 1.5, 3.0)),
        (atanh, (0.0, 0.3, -0.4)))

    @testset "TaylorN and Taylor1{Taylor1} (2 variables)" begin
        sp = JetSpace(order = 10, variables = [:τ, :σ])
        x, y = variables(sp)
        set_taylor1_varname(2, ["t", "s"])
        t = Taylor1(5)                          # inner variable (τ)
        s = Taylor1([zero(t), one(t)], 5)       # outer variable (σ)
        for (f, x0s) in _cmp_funs()
            @testset "$f" begin
                for x0 in x0s
                    @test _same_coeffs(f(x0 + x + y), f(x0 + t + s))
                    @test _same_coeffs(f(x0 + 0.5x - 0.3y), f(x0 + 0.5t - 0.3s))
                end
            end
        end
    end

    @testset "TaylorN and Taylor1{Taylor1{Taylor1}} (3 variables)" begin
        sp = JetSpace(order = 9, variables = [:τ, :σ, :ρ])
        x, y, z = variables(sp)
        set_taylor1_varname(3, ["t", "s", "r"])
        t = Taylor1(3)                          # innermost variable (τ)
        s = Taylor1([zero(t), one(t)], 3)       # middle variable (σ)
        r = Taylor1([zero(s), one(s)], 3)       # outermost variable (ρ)
        for (f, x0s) in _cmp_funs()
            @testset "$f" begin
                for x0 in x0s
                    @test _same_coeffs(f(x0 + x + y + z), f(x0 + t + s + r))
                    @test _same_coeffs(f(x0 + 0.5x - 0.3y + 0.2z), f(x0 + 0.5t - 0.3s + 0.2r))
                end
            end
        end
    end

    @testset "TaylorN and Taylor1{TaylorN} (2 combined vars each)" begin
        sp = JetSpace(order = 10, variables = [:τ, :σ])
        x, y = variables(sp)
        sp1 = JetSpace(order = 5, variables = [:τ])
        τ = variables(sp1)[1]                   # inner level: TaylorN in one variable (τ)
        s = Taylor1([zero(τ), one(τ)], 5)       # outer level: Taylor1 variable (σ)
        # mixed arithmetic: a TaylorN behaves as the constant Taylor1{TaylorN} it represents
        b = 1.0 + 0.5τ
        bT = Taylor1(b, 5)
        for op in (+, -, *, /)
            @test op(b, 2 + s) == op(bT, 2 + s)
            @test op(2 + s, b) == op(2 + s, bT)
        end
        for (f, x0s) in _cmp_funs()
            @testset "$f" begin
                for x0 in x0s
                    @test _same_coeffs(f(x0 + x + y), f(x0 + τ + s))
                    @test _same_coeffs(f(x0 + 0.5x - 0.3y), f(x0 + 0.5τ - 0.3s))
                end
            end
        end
    end

    @testset "TaylorN Taylor1{TaylorN} (3 combined vars each)" begin
        sp3 = JetSpace(order = 8, variables = [:τ, :σ, :ρ])
        x, y, z = variables(sp3)
        sp2 = JetSpace(order = 4, variables = [:τ, :σ])
        u, v = variables(sp2)                   # inner level: TaylorN in two variables (τ, σ)
        r = Taylor1([zero(u), one(u)], 4)       # outer level: Taylor1 variable (ρ)
        for (f, x0s) in _cmp_funs()
            @testset "$f" begin
                for x0 in x0s
                    @test _same_coeffs(f(x0 + x + y + z), f(x0 + u + v + r))
                    # (subtracting the outer variable: TaylorN - Taylor1{TaylorN})
                    @test _same_coeffs(f(x0 + 0.5x + 0.2y - 0.3z), f(x0 + 0.5u + 0.2v - 0.3r))
                end
            end
        end
    end
end

# Integrity checks of the expansions of the elementary functions using: (i) compare the
# evaluation at a small increment of the expansion around x0, f(x0 + h), and
# computed by Base; and (ii) the chain rule: the derivative of f(x) equals f′(x) x′
@testset "Evaluation and chain rule tests for the elementary functions" begin
    # function, its derivative, and three expansion points within its domain
    funs = ((exp, exp, (0.0, 0.3, -0.4)), (expm1, exp, (0.0, 0.3, -0.4)),
        (log, x -> 1 / x, (1.0, 1.5, 0.7)), (log1p, x -> 1 / (1 + x), (0.0, 0.3, -0.4)),
        (sqrt, x -> 1 / (2 * sqrt(x)), (1.0, 1.5, 0.7)), (sin, cos, (0.0, 0.3, -0.4)),
        (cos, x -> -sin(x), (0.0, 0.3, -0.4)), (tan, x -> 1 + tan(x)^2, (0.0, 0.3, -0.4)),
        (sinpi, x -> pi * cospi(x), (0.0, 0.3, -0.4)),
        (cospi, x -> -pi * sinpi(x), (0.0, 0.3, -0.4)),
        (asin, x -> 1 / sqrt(1 - x^2), (0.0, 0.3, -0.4)),
        (acos, x -> -1 / sqrt(1 - x^2), (0.0, 0.3, -0.4)),
        (atan, x -> 1 / (1 + x^2), (0.0, 0.3, -0.4)), (sinh, cosh, (0.0, 0.3, -0.4)),
        (cosh, sinh, (0.0, 0.3, -0.4)), (tanh, x -> 1 - tanh(x)^2, (0.0, 0.3, -0.4)),
        (asinh, x -> 1 / sqrt(x^2 + 1), (0.0, 0.3, -0.4)),
        (acosh, x -> 1 / sqrt(x^2 - 1), (2.0, 1.5, 3.0)),
        (atanh, x -> 1 / (1 - x^2), (0.0, 0.3, -0.4)))

    @testset "Evaluation at a small increment" begin
        sp = JetSpace(order = 12, variables = [:τ, :σ])
        x, y = variables(sp)
        t = Taylor1(12)
        s = Taylor1([zero(t), one(t)], 12)          # Taylor1{Taylor1}: σ outer, τ inner
        sp1 = JetSpace(order = 12, variables = [:τ])
        τ = variables(sp1)[1]
        sτ = Taylor1([zero(τ), one(τ)], 12)         # Taylor1{TaylorN}: σ outer, τ inner
        h1, h2 = 0.011, -0.007
        for (f, _, x0s) in funs, x0 in x0s
            exact = f(x0 + h1 + h2)
            @test isapprox(evaluate(f(x0 + t), h1 + h2), exact; rtol = 32 * eps())
            @test isapprox(evaluate(f(x0 + x + y), [h1, h2]), exact; rtol = 32 * eps())
            @test isapprox(evaluate(evaluate(f(x0 + t + s), h2), h1), exact; rtol = 32 * eps())
            @test isapprox(evaluate(evaluate(f(x0 + τ + sτ), h2), [h1]), exact; rtol = 32 * eps())
        end
    end

    @testset "Chain rule" begin
        # coefficients to compare: all for Taylor1; those of degree < order for TaylorN
        # (`differentiate` keeps the order of a TaylorN, and f′(x) x′ is incomplete at that degree)
        coeffs_cmp(p::Taylor1) = collect(p.coeffs)
        coeffs_cmp(p::TaylorN) =
            reduce(vcat, (collect(p.coeffs[d+1].coeffs) for d in 0:TS.order(p)-1))
        # the same series up to rounding, relative to its largest coefficient (a small
        # coefficient may carry the rounding error of large terms that cancel)
        function same_series(a, b; rtol = 32 * eps())
            ca, cb = coeffs_cmp(a), coeffs_cmp(b)
            scale = max(maximum(abs, ca), maximum(abs, cb))
            return length(ca) == length(cb) &&
                all(abs(p - q) ≤ rtol * scale for (p, q) in zip(ca, cb))
        end
        sp = JetSpace(order = 10, variables = [:u, :v])
        u, v = variables(sp)
        t = Taylor1(15)
        for (f, df, x0s) in funs, x0 in x0s
            xt = x0 + t + 0.5 * t^2
            @test same_series(differentiate(f(xt)), df(xt) * differentiate(xt))
            xN = x0 + u - 0.3v + 0.5u * v
            @test same_series(differentiate(f(xN), 1), df(xN) * differentiate(xN, 1))
            @test same_series(differentiate(f(xN), 2), df(xN) * differentiate(xN, 2))
        end
    end
end