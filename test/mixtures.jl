# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries

using Test
# using LinearAlgebra

@testset "Tests with mixtures of Taylor1 and TaylorN" begin
    set_taylor1_varname(1, "t")
    @test TS.NumberNotSeries == Union{Real,Complex}
    @test TS.NumberNotSeriesN == Union{Real,Complex,Taylor1}

    set_variables("x", numvars=2, order=6)
    xH = HomogeneousPolynomial(Int, 1)
    yH = HomogeneousPolynomial(Int, 2)
    tN = Taylor1(TaylorN{Float64}, 3)
    @test findfirst(tN) == 1

    @test convert(eltype(tN), tN) == tN
    @test eltype(xH) == HomogeneousPolynomial{Int}
    @test TS.numtype(xH) == Int
    @test eltype(tN) == Taylor1{TaylorN{Float64}}
    @test TS.numtype(tN) == TaylorN{Float64}
    @test normalize_taylor(tN) == tN
    @test tN.order == 3

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
    @test get_order(HomogeneousPolynomial([Taylor1(1), 1.0+Taylor1(2)])) == 1
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
        @test tN1.coeffs[i].order == zeroN1.coeffs[i].order == oneN1.coeffs[i].order
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
    v = zeros(TaylorN{Float64},2)
    @test isnothing(evaluate!([t1N, t1N^2], 0.0, v))
    @test v == [TaylorN(1), TaylorN(1)^2]
    @test tN1() == t
    @test evaluate(tN1, :x₁ => 1.0) == TaylorN([HomogeneousPolynomial([1.0+t]), zero(xHt), yHt^2])
    @test evaluate(tN1, 1, 1.0) == TaylorN([HomogeneousPolynomial([1.0+t]), zero(xHt), yHt^2])
    @test evaluate(t, t1N) == t1N
    @test evaluate(t1N, 0.5) == t1N[0] + t1N[1]/2 + t1N[2]/4
    @test evaluate(t1N, [t1N[0], zero(t1N[0])]) == Taylor1([t1N[0], t1N[1]], t.order)
    @test evaluate(t1N, 2, 0.0) == Taylor1([t1N[0], t1N[1]], t.order)
    @test evaluate(t1N, 1, 0.0) == Taylor1([zero(t1N[0]), t1N[1], t1N[2]], t.order)

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
    ee = Taylor1(t1N[0:5], 6)
    for ord in eachindex(t1N)
        TS.differentiate!(ee, exp(t1N), ord)
    end
    @test iszero(ee[6])
    @test getcoeff.(ee, 0:5) == getcoeff.(exp(t1N), 0:5)
    ee = differentiate(t1N, get_order(t1N))
    @test iszero(ee)
    @test iszero(get_order(ee))

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
    xN, yN = get_variables()
    @test evaluate(1.0 + t + t^2, xN) == 1.0 + xN + xN^2
    v1 = [1.0 + t + t^2 + t^4, 1.0 - t^2 + t^3]
    @test v1(yN^2) == [1.0 + yN^2 + yN^4, 1.0 - yN^4 + yN^6]
    tN = Taylor1([zero(xN), one(xN)], 4)
    q1N = 1 + yN*tN + xN*tN^4
    @test q1N(-1.0) == 1.0 - yN + xN
    @test q1N(-xN^2) == 1.0 - xN^2*yN

    # See #92 and #94
    δx, δy = set_variables("δx δy")
    xx = 1+Taylor1(δx, 5)
    yy = 1+Taylor1(δy, 5)
    tt = Taylor1([zero(δx), one(δx)], xx.order)
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
    @test sqrt(xx) == Taylor1(sqrt(xx[0]), xx.order)
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
    dx = set_variables("x", numvars=4, order=10)
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

    X, Y = set_variables(Taylor1{Float64}, "x y")
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

    P = rndTN(get_order(), 3)
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

    X, Y = set_variables(BigFloat, "x y", numvars=2, order=6)
    p1N = Taylor1([X^2,X*Y,Y+X,Y^2])
    q1N = Taylor1([X^2,(1.0+sqrt(eps(BigFloat)))*X*Y,Y+X,Y^2])
    @test p1N ≈ p1N
    @test p1N ≈ q1N

    Pv = [rndTN(get_order(), 3), rndTN(get_order(), 3)]
    Qv = convert.(Taylor1{TaylorN{Float64}}, Pv)

    @test TS.jacobian(Pv) == TS.jacobian(Qv)

    @test_throws ArgumentError Taylor1(2) + TaylorN(1)
    @test_throws ArgumentError Taylor1(2) - TaylorN(1)
    @test_throws ArgumentError Taylor1(2) * TaylorN(1)
    @test_throws ArgumentError TaylorN(2) / Taylor1(1)

    # Issue #342 and PR #343
    z0N = -1.333+get_variables()[1]
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

    a = sum(exp.(get_variables()).^2)
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

    dq = get_variables()
    x = Taylor1(exp.(dq), 5)
    x[1] = sin(dq[1]*dq[2])
    @test x[1] == sin(dq[1]*dq[2])
    @test x[1] !== sin(dq[1]*dq[2])

    @testset "Test Base.float overloads for Taylor1 and TaylorN mixtures" begin
        q = get_variables(Int)
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
    @test titii / tii == ti
    @test get_order(titii/tii) == get_order(tii)-1
    @test titii / ti == tii
    @test get_order(titii/ti) == get_order(tii)
    @test ti^2-tii^2 == (ti+tii)*(ti-tii)
    @test findfirst(ti^2-tii^2) == 0
    @test sin(tii) ≈ Taylor1(one(ti) .* sin(Taylor1(10)).coeffs, 9)
    @test tii(1 + ti) == 1 + ti
    @test tii(1 + ti) isa Taylor1{Float64}
    @test ti(1 + tii) == 1 + tii
    @test constant_term(ti+tii) == ti
    @test linear_polynomial(ti*tii) == Taylor1([zero(ti), ti], 9)
    @test get_order(linear_polynomial(tii)) == get_order(tii)
    @test nonlinear_polynomial(tii+ti*tii^2) == Taylor1([zero(ti), zero(ti), ti], 9)
    @test ti(1 + tii) isa Taylor1{Taylor1{Float64}}
    @test sqrt(titii^2) == titii
    @test get_order(sqrt(titii^2)) == get_order(tii) >> 1
    @test (titii^3)^(1/3) == titii
    @test get_order(sqrt(titii^2)) == get_order(tii) >> 1
    ti2tii = ti^2 * tii
    tti = (ti2tii/tii)/ti
    @test get_order(tti) == get_order(tii)-1
    @test get_order(tti[0]) == get_order(ti)-1
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

    # Back to default
    set_taylor1_varname(1, "t")
    @test TS._params_Taylor1_.var_name == ["t"]
end
