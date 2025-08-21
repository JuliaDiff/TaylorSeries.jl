# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries, IntervalArithmetic

using Test
# eeuler = Base.MathConstants.e

setdisplay(:full)

@testset "Tests Taylor1 and TaylorN expansions over Intervals" begin
    a = interval(1, 2)
    b = interval(-1, 1)
    c = interval(0, 1)
    ithree = interval(3)
    ifour = interval(4)
    ifive = interval(5)
    isix = interval(6)
    iten = interval(10)
    p3(x, a) = x^3 + ithree*a*x^2 + ithree*a^2*x + a^3
    p4(x, a) = x^4 + ifour*a*x^3 + isix*a^2*x^2 + ifour*a^3*x + a^4
    p5(x, a) = x^5 + ifive*a*x^4 + iten*a^2*x^3 + iten*a^3*x^2 + ifive*a^4*x + a^5

    ti = Taylor1(Interval{Float64}, 10)
    x, y = set_variables(Interval{Float64}, "x y")

    @test eltype(ti) == Taylor1{Interval{Float64}}
    @test eltype(x) == TaylorN{Interval{Float64}}
    @test TS.numtype(ti) == Interval{Float64}
    @test TS.numtype(x) == Interval{Float64}
    @test normalize_taylor(ti) == ti
    @test normalize_taylor(x) == x

    @test ti == Taylor1(typeof(c), 3)
    @test p3(ti,-a) == (ti-a)^3 == (ti-a)^3.0
    @test p4(ti,-a) == (ti-a)^4 == (ti-a)^4.0
    @test p5(ti,-a) == (ti-a)^5 == (ti-a)^5.0
    @test (ti-b)^3 == (ti-b)^3.0
    @test all(issubset_interval.((p3(ti,-b)).coeffs, ((ti-b)^3).coeffs))
    @test p4(ti,-b) == (ti-b)^4 == (ti-b)^4.0
    @test (ti-b)^5 == (ti-b)^5.0
    @test all(issubset_interval.((p5(ti,-b)).coeffs, ((ti-b)^5).coeffs))
    @test ti^(4//1) == ti^4

    @test x[1] !== y[0]
    @test p3(x,-y) == (x-y)^3 == (x-y)^3.0
    @test p4(x,-y) == (x-y)^4 == (x-y)^4.0
    @test p5(x,-y) == (x-y)^5 == (x-y)^5.0
    @test p3(x,-a) == (x-a)^3 == (x-a)^3.0
    @test p4(x,-a) == (x-a)^4 == (x-a)^4.0
    @test p5(x,-a) == (x-a)^5 == (x-a)^5.0
    @test (x-b)^3 == (x-b)^3.0
    for ind in eachindex(p3(x,-b))
        @test all(issubset_interval.((p3(x,-b)[ind]).coeffs, (((x-b)^3)[ind]).coeffs))
    end
    @test p4(x,-b) == (x-b)^4
    @test (x-b)^5 == (x-b)^5.0
    r1 = p5(x,-b)
    r2 = (x-b)^5
    for ind in eachindex(p5(x,-b))
        @test all(issubset_interval.(r1[ind].coeffs, r2[ind].coeffs))
        @test all(isguaranteed.(getfield(r1[ind], :coeffs)))
        @test all(isguaranteed.(getfield(r2[ind], :coeffs)))
    end
    @test y^(4//1) == y^4


    # Tests `evaluate`
    @test isequal_interval(evaluate(p3(x,y), [a, -b]), p3(a, -b))
    @test isequal_interval(evaluate(p4(x,y), [a, -b]), p4(a, -b))
    @test isequal_interval((p5(x,y))([a, b]), p5(a, b))
    @test issubset_interval((a-b)^4, ((x-y)^4)([a, b]))
    @test isequal_interval((((x-y)^4)[4])([a, b]), interval(-39, 81))

    p4n = normalize_taylor(p4(x,y), [a, b], true)
    @test issubset_interval(interval(0, 16), p4n([b, b]))
    p5n = normalize_taylor(p5(x,y), [a, b], true)
    @test issubset_interval(interval(-32, 32), p5n([b, b]))

    p4n = normalize_taylor(p4(x,y), [a, b], false)
    @test issubset_interval(interval(0, 16), p4n([c, c]))
    p5n = normalize_taylor(p5(x,y), [a, b], false)
    @test issubset_interval(interval(0, 32), p5n([c, c]))

    @test isequal_interval(evaluate(x*y^3, [b, b]), b)
    @test isequal_interval(evaluate(x*y^2, [b, b]), b)
    @test isequal_interval(evaluate(x^2*y^2, [b, b]), c)

    ii = b
    t = Taylor1(1)
    @test issubset_interval(interval(0, 2), (1+t)(ii))
    t = Taylor1(2)
    @test issubset_interval(interval(0, 4), ((1+t)^2)(ii))

    ii = interval(0, 6)
    t = Taylor1(4)
    f(x) = 0.1 * x^3 - 0.5 * x^2 + 1
    ft = f(t)
    f1 = normalize_taylor(ft, ii, true)
    f2 = normalize_taylor(ft, ii, false)
    @test !isguaranteed(f(ii))
    @test !isguaranteed(ft(ii))
    @test !isguaranteed(f1(b))
    @test !isguaranteed(f2(c))
    @test issubset_interval(interval(-23/27, f(6)), f(ii))
    @test issubset_interval(interval(-23/27, f(6)), ft(ii))
    @test issubset_interval(interval(-23/27, f(6)), f1(b))
    @test issubset_interval(interval(-23/27, f(6)), f2(c))
    @test issubset_interval(f1(b), f(ii))
    @test diam(f1(b)) < diam(f2(c))


    # An example from Makino's thesis
    ii = c
    t = Taylor1(5)
    g(x) = 1 - x^4 + x^5
    gt1 = g(t)
    gn1 = normalize_taylor(gt1, c, true)
    @test issubset_interval(interval(g(4/5), 1), g(ii))
    @test issubset_interval(interval(g(4/5), 1), gt1(ii))
    @test issubset_interval(interval(g(4/5), 1), gn1(b))
    @test isinterior(gn1(b), g(ii))
    @test diam(gn1(b)) < diam(gt1(ii))
    #
    ti = Taylor1(typeof(c), 5)
    gg(x) = interval(1) - x^4 + x^5
    gt = gg(ti)
    @test all(isguaranteed.(getfield(gt, :coeffs)))
    gn2 = normalize_taylor(gt, c, true)
    @test all(isguaranteed.(getfield(gn2, :coeffs)))
    @test issubset_interval(interval(g(4/5), 1), gg(ii))
    @test issubset_interval(interval(g(4/5), 1), gt(ii))
    @test issubset_interval(interval(g(4/5), 1), gn2(b))
    @test isinterior(gn2(b), g(ii))
    @test diam(gn2(b)) < diam(gt(ii))


    # Test display for Taylor1{Complex{Interval{T}}}
    vc = [complex(interval(1.5, 2), interval(0, 0)),
        complex(interval(-2, -1), interval(-1, 1 )),
        complex(interval(-1, 1.5), interval(-1, 1.5)),
        complex(interval(0, 0), interval(-1, 1.5))]
    displayBigO(false)
    @test string(Taylor1(vc, 5)) ==
        " ( Interval{Float64}(1.5, 2.0, com, true) + im*Interval{Float64}(0.0, 0.0, com, true) ) +" *
        " ( Interval{Float64}(-2.0, -1.0, com, true) + im*Interval{Float64}(-1.0, 1.0, com, true) ) t" *
        " + ( Interval{Float64}(-1.0, 1.5, com, true) + im*Interval{Float64}(-1.0, 1.5, com, true) ) t²" *
        " + ( Interval{Float64}(0.0, 0.0, com, true) + im*Interval{Float64}(-1.0, 1.5, com, true) ) t³ "
    displayBigO(true)
    @test string(Taylor1(vc, 5)) ==
        " ( Interval{Float64}(1.5, 2.0, com, true) + im*Interval{Float64}(0.0, 0.0, com, true) ) +" *
        " ( Interval{Float64}(-2.0, -1.0, com, true) + im*Interval{Float64}(-1.0, 1.0, com, true) ) t" *
        " + ( Interval{Float64}(-1.0, 1.5, com, true) + im*Interval{Float64}(-1.0, 1.5, com, true) ) t²" *
        " + ( Interval{Float64}(0.0, 0.0, com, true) + im*Interval{Float64}(-1.0, 1.5, com, true) ) t³ + 𝒪(t⁶)"


    # Iss 351 (inspired by a test in ReachabilityAnalysis)
    p1 = Taylor1([interval(0, 0), interval(0, 0.1) + interval(0, 0.01) * y], 4)
    p2 = Taylor1([interval(0, 0), interval(0, 0.5) + interval(0, 0.02) * x + interval(0, 0.03) * y], 4)
    @test evaluate([p1, p2], interval(0, 1)) == [p1[1], p2[1]]
    @test typeof(p1(interval(0, 1))) == TaylorN{Interval{Float64}}

    # Tests related to Iss #311
    # `sqrt` and `pow` defined on Interval(0,Inf)
    @test_throws DomainError sqrt(ti)
    @test_throws DomainError ti^(1/4)
    @test sqrt(interval(0.0, 1.e-15) + ti) == sqrt(interval(-1.e-15, 1.e-15) + ti)
    aa = sqrt(sqrt(interval(0.0, 1.e-15) + ti))
    @test aa^0 == one(aa)
    @test aa^1 == aa
    @test !(aa^1 === aa)
    @test aa == sqrt(sqrt(interval(-1.e-15, 1.e-15) + ti))
    bb = (interval(0.0, 1.e-15) + ti)^(1/4)
    @test bb == (interval(-1.e-15, 1.e-15) + ti)^(1/4)
    @test all(isinterior.(aa.coeffs[2:end], bb.coeffs[2:end]))
    @test (one(ti)+ti)^-1 == one(ti) - ti + ti^2 - ti^3 + ti^4 - ti^5
    @test_throws DomainError sqrt(x)
    @test_throws DomainError x^(1/4)
    @test sqrt(interval(-1,1)+x) == sqrt(interval(0,1)+x)
    @test (interval(-1,1)+x)^(1/4) == (interval(0,1)+x)^(1/4)
    @test (one(y)+y)^-1 == one(y) - y + y^2 - y^3 + y^4 - y^5 + y^6

    #
    @test all(issubset_interval.(one(ti)[:], (sin(ti)^2+cos(ti)^2)[:]))
    @test tan(ti)^2 == sec(ti)^2 - one(ti)
    for ind in eachindex(one(x))
        @test all(issubset_interval.((one(x)[ind]).coeffs, (sin(x)^2+cos(x)^2)[ind].coeffs))
    end

    # `log` defined on Interval(0,Inf)
    @test_throws DomainError log(ti)
    @test log(interval(0.0, 1.e-15) + ti) == log(interval(-1.e-15, 1.e-15) + ti)
    @test_throws DomainError log(y)
    @test log(interval(0.0, 1.e-15) + y) == log(interval(-1.e-15, 1.e-15) + y)
    @test all(in_interval.(log(0.875+Taylor1(5))[:], log1p(-0.125+ti)[:]))

    # `asin` and `acos` defined on interval(-1,1)
    @test_throws DomainError asin(interval(1.0, 2.0) + ti)
    @test asin(interval(-2.0, 0.0) + ti) == asin(interval(-1,0) + ti)
    @test_throws DomainError acos(interval(1.0, 2.0) + ti)
    @test acos(interval(-2.0, 0.0) + ti) == acos(interval(-1,0) + ti)
    @test_throws DomainError asin(interval(1.0, 2.0) + x)
    @test asin(interval(-2.0, 0.0) + x) == asin(interval(-1,0) + x)
    @test_throws DomainError acos(interval(1.0, 2.0) + x)
    @test acos(interval(-2.0, 0.0) + x) == acos(interval(-1,0) + x)
    # acosh defined on interval(1,Inf)
    @test_throws DomainError acosh(interval(0.0, 1.0) + ti)
    @test acosh(interval(0.0, 2.0) + ti) == acosh(interval(1.0, 2.0) + ti)
    @test_throws DomainError acosh(interval(0.0, 1.0) + x)
    @test acosh(interval(0.0, 2.0) + x) == acosh(interval(1.0, 2.0) + x)
    # atanh defined on interval(-1,1)
    @test_logs (:warn, "ill-formed bare interval [a, b] with a = Inf, b = Inf. Empty interval is returned") @test_throws DomainError atanh(interval(1.0, 1.0) + ti)
    @test atanh(interval(-2.0, 0.0) + ti) == atanh(interval(-1.0, 0.0) + ti)
    @test_logs (:warn, "ill-formed bare interval [a, b] with a = Inf, b = Inf. Empty interval is returned") @test_throws DomainError atanh(interval(1.0, 1.0) + y)
    @test atanh(interval(-2.0, 0.0) + y) == atanh(interval(-1.0, 0.0) + y)
end
