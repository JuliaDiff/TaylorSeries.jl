# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries, IntervalArithmetic

using Test
eeuler = Base.MathConstants.e

@testset "Tests Taylor1 and TaylorN expansions over Intervals" begin
    a = 1..2
    b = -1 .. 1
    p4(x, a) = x^4 + 4*a*x^3 + 6*a^2*x^2 + 4*a^3*x + a^4
    p5(x, a) = x^5 + 5*a*x^4 + 10*a^2*x^3 + 10*a^3*x^2 + 5*a^4*x + a^5

    ti = Taylor1(Interval{Float64}, 10)
    x, y = set_variables(Interval{Float64}, "x y")

    @test eltype(ti) == Interval{Float64}
    @test eltype(x) == Interval{Float64}

    @test p4(ti,-a) == (ti-a)^4
    @test p5(ti,-a) == (ti-a)^5
    @test p4(ti,-b) == (ti-b)^4
    @test all((p5(ti,-b)).coeffs .⊆ ((ti-b)^5).coeffs)


    @test p4(x,-y) == (x-y)^4
    @test p5(x,-y) == (x-y)^5
    @test p4(x,-a) == (x-a)^4
    @test p5(x,-a) == (x-a)^5
    @test p4(x,-b) == (x-b)^4
    for ind in eachindex(p5(x,-b))
        @test all((p5(x,-b)[ind]).coeffs .⊆ (((x-b)^5)[ind]).coeffs)
    end

    # Tests `evaluate`
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
end
