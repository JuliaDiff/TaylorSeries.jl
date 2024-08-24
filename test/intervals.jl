# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries, IntervalArithmetic

using Test
# eeuler = Base.MathConstants.e

@testset "Tests Taylor1 and TaylorN expansions over Intervals" begin
    a = 1..2
    b = -1 .. 1
    p3(x, a) = x^3 + 3*a*x^2 + 3*a^2*x + a^3
    p4(x, a) = x^4 + 4*a*x^3 + 6*a^2*x^2 + 4*a^3*x + a^4
    p5(x, a) = x^5 + 5*a*x^4 + 10*a^2*x^3 + 10*a^3*x^2 + 5*a^4*x + a^5

    ti = Taylor1(Interval{Float64}, 10)
    x, y = set_variables(Interval{Float64}, "x y")

    # @test eltype(ti) == Interval{Float64}
    # @test eltype(x) == Interval{Float64}
    @test eltype(ti) == Taylor1{Interval{Float64}}
    @test eltype(x) == TaylorN{Interval{Float64}}
    @test TS.numtype(ti) == Interval{Float64}
    @test TS.numtype(x) == Interval{Float64}
    @test normalize_taylor(ti) == ti
    @test normalize_taylor(x) == x

    @test p3(ti,-a) == (ti-a)^3 == (ti-a)^3.0
    @test p4(ti,-a) == (ti-a)^4 == (ti-a)^4.0
    @test p5(ti,-a) == (ti-a)^5 == (ti-a)^5.0
    @test (ti-b)^3 == (ti-b)^3.0
    @test all((p3(ti,-b)).coeffs .⊆ ((ti-b)^3).coeffs)
    @test p4(ti,-b) == (ti-b)^4 == (ti-b)^4.0
    @test (ti-b)^5 == (ti-b)^5.0
    @test all((p5(ti,-b)).coeffs .⊆ ((ti-b)^5).coeffs)


    @test p3(x,-y) == (x-y)^3 == (x-y)^3.0
    @test p4(x,-y) == (x-y)^4 == (x-y)^4.0
    @test p5(x,-y) == (x-y)^5 == (x-y)^5.0
    @test p3(x,-a) == (x-a)^3 == (x-a)^3.0
    @test p4(x,-a) == (x-a)^4 == (x-a)^4.0
    @test p5(x,-a) == (x-a)^5 == (x-a)^5.0
    @test (x-b)^3 == (x-b)^3.0
    for ind in eachindex(p3(x,-b))
        @test all((p3(x,-b)[ind]).coeffs .⊆ (((x-b)^3)[ind]).coeffs)
    end
    @test p4(x,-b) == (x-b)^4
    @test (x-b)^5 == (x-b)^5.0
    for ind in eachindex(p5(x,-b))
        @test all((p5(x,-b)[ind]).coeffs .⊆ (((x-b)^5)[ind]).coeffs)
    end

    # Tests `evaluate`
    @test evaluate(p3(x,y), IntervalBox(a,-b)) == p3(a, -b)
    @test evaluate(p4(x,y), IntervalBox(a,-b)) == p4(a, -b)
    @test (p5(x,y))(IntervalBox(a,b)) == p5(a, b)
    @test (a-b)^4 ⊆ ((x-y)^4)(a × b)
    @test (((x-y)^4)[4])(a × b) == -39 .. 81

    p4n = normalize_taylor(p4(x,y), a × b, true)
    @test (0..16) ⊆ p4n((-1..1)×(-1..1))
    p5n = normalize_taylor(p5(x,y), a × b, true)
    @test (-32 .. 32) ⊆ p5n((-1..1)×(-1..1))

    p4n = normalize_taylor(p4(x,y), a × b, false)
    @test (0..16) ⊆ p4n((0..1)×(0..1))
    p5n = normalize_taylor(p5(x,y), a × b, false)
    @test (0..32) ⊆ p5n((0..1)×(0..1))

    @test evaluate(x*y^3, (-1..1)×(-1..1)) == (-1..1)
    @test evaluate(x*y^2, (-1..1)×(-1..1)) == (-1..1)
    @test evaluate(x^2*y^2, (-1..1)×(-1..1)) == (0..1)

    ii = -1..1
    t = Taylor1(1)
    @test 0..2 ⊆ (1+t)(ii)
    t = Taylor1(2)
    @test 0..4 ⊆ ((1+t)^2)(ii)

    ii = 0..6
    t = Taylor1(4)
    f(x) = 0.1 * x^3 - 0.5*x^2 + 1
    ft = f(t)
    f1 = normalize_taylor(ft, ii, true)
    f2 = normalize_taylor(ft, ii, false)
    @test Interval(-23/27, f(6)) ⊆ f(ii)
    @test Interval(-23/27, f(6)) ⊆ ft(ii)
    @test Interval(-23/27, f(6)) ⊆ f1(-1..1)
    @test Interval(-23/27, f(6)) ⊆ f2(0..1)
    @test f1(-1..1) ⊆ f(ii)
    @test diam(f1(-1..1)) < diam(f2(0..1))

    # An example from Makino's thesis
    ii = 0..1
    t = Taylor1(5)
    g(x) = 1 - x^4 + x^5
    gt = g(t)
    g1 = normalize_taylor(gt, 0..1, true)
    @test Interval(g(4/5),1) ⊆ g(ii)
    @test Interval(g(4/5),1) ⊆ gt(ii)
    @test Interval(g(4/5),1) ⊆ g1(-1..1)
    @test g1(-1..1) ⊂ g(ii)
    @test diam(g1(-1..1)) < diam(gt(ii))

    # Test display for Taylor1{Complex{Interval{T}}}
    vc = [complex(1.5 .. 2, 0..0 ), complex(-2  .. -1, -1 .. 1 ),
        complex( -1 .. 1.5, -1 .. 1.5), complex( 0..0, -1 .. 1.5)]
    displayBigO(false)
    @test string(Taylor1(vc, 5)) ==
        " ( [1.5, 2] + [0, 0]im ) - ( [1, 2] + [-1, 1]im ) t + ( [-1, 1.5] + [-1, 1.5]im ) t² + ( [0, 0] + [-1, 1.5]im ) t³ "
    displayBigO(true)
    @test string(Taylor1(vc, 5)) ==
        " ( [1.5, 2] + [0, 0]im ) - ( [1, 2] + [-1, 1]im ) t + ( [-1, 1.5] + [-1, 1.5]im ) t² + ( [0, 0] + [-1, 1.5]im ) t³ + 𝒪(t⁶)"

    # Iss 351 (inspired by a test in ReachabilityAnalysis)
    p1 = Taylor1([0 .. 0, (0 .. 0.1) + (0 .. 0.01) * y], 4)
    p2 = Taylor1([0 .. 0, (0 .. 0.5) + (0 .. 0.02) * x + (0 .. 0.03) * y], 4)
    @test evaluate([p1, p2], 0 .. 1) == [p1[1], p2[1]]
    @test typeof(p1(0 .. 1)) == TaylorN{Interval{Float64}}

    # Tests related to Iss #311
    # `sqrt` and `pow` defined on Interval(0,Inf)
    @test_throws DomainError sqrt(ti)
    @test sqrt(Interval(0.0, 1.e-15) + ti) == sqrt(Interval(-1.e-15, 1.e-15) + ti)
    aa = sqrt(sqrt(Interval(0.0, 1.e-15) + ti))
    @test aa == sqrt(sqrt(Interval(-1.e-15, 1.e-15) + ti))
    bb = (Interval(0.0, 1.e-15) + ti)^(1/4)
    @test bb == (Interval(-1.e-15, 1.e-15) + ti)^(1/4)
    @test all(aa.coeffs[2:end] .⊂ bb.coeffs[2:end])
    @test_throws DomainError sqrt(x)
    @test sqrt(Interval(-1,1)+x) == sqrt(Interval(0,1)+x)
    @test (Interval(-1,1)+x)^(1/4) == (Interval(0,1)+x)^(1/4)

    # `log` defined on Interval(0,Inf)
    @test_throws DomainError log(ti)
    @test log(Interval(0.0, 1.e-15) + ti) == log(Interval(-1.e-15, 1.e-15) + ti)
    @test_throws DomainError log(y)
    @test log(Interval(0.0, 1.e-15) + y) == log(Interval(-1.e-15, 1.e-15) + y)
    # `asin` and `acos` defined on Interval(-1,1)
    @test_throws DomainError asin(Interval(1.0 .. 2.0) + ti)
    @test asin(Interval(-2.0 .. 0.0) + ti) == asin(Interval(-1,0) + ti)
    @test_throws DomainError acos(Interval(1.0 .. 2.0) + ti)
    @test acos(Interval(-2.0 .. 0.0) + ti) == acos(Interval(-1,0) + ti)
    @test_throws DomainError asin(Interval(1.0 .. 2.0) + x)
    @test asin(Interval(-2.0 .. 0.0) + x) == asin(Interval(-1,0) + x)
    @test_throws DomainError acos(Interval(1.0 .. 2.0) + x)
    @test acos(Interval(-2.0 .. 0.0) + x) == acos(Interval(-1,0) + x)
    # acosh defined on Interval(1,Inf)
    @test_throws DomainError acosh(Interval(0.0 .. 1.0) + ti)
    @test acosh(Interval(0.0 .. 2.0) + ti) == acosh(Interval(1.0 .. 2.0) + ti)
    @test_throws DomainError acosh(Interval(0.0 .. 1.0) + x)
    @test acosh(Interval(0.0 .. 2.0) + x) == acosh(Interval(1.0 .. 2.0) + x)
    # atanh defined on Interval(-1,1)
    @test_throws DomainError atanh(Interval(1.0 .. 1.0) + ti)
    @test atanh(Interval(-2.0 .. 0.0) + ti) == atanh(Interval(-1.0 .. 0.0) + ti)
    @test_throws DomainError atanh(Interval(1.0 .. 1.0) + y)
    @test atanh(Interval(-2.0 .. 0.0) + y) == atanh(Interval(-1.0 .. 0.0) + y)
end
