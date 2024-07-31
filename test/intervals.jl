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
    # p4(x, a) = x^4 + 4*a*x^3 + 6*a^2*x^2 + 4*a^3*x + a^4
    # p5(x, a) = x^5 + 5*a*x^4 + 10*a^2*x^3 + 10*a^3*x^2 + 5*a^4*x + a^5
    ifour = interval(4)
    ifive = interval(5)
    isix = interval(6)
    iten = interval(10)
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

    @test p4(ti,-a) == (ti-a)^4
    @test p5(ti,-a) == (ti-a)^5
    @test p4(ti,-b) == (ti-b)^4
    @test all(issubset_interval.((p5(ti,-b)).coeffs, ((ti-b)^5).coeffs))

    @test p4(x,-y) == (x-y)^4
    @test p5(x,-y) == (x-y)^5
    @test p4(x,-a) == (x-a)^4
    @test p5(x,-a) == (x-a)^5
    @test p4(x,-b) == (x-b)^4
    r1 = p5(x,-b)
    r2 = (x-b)^5
    for ind in eachindex(p5(x,-b))
        @test all(issubset_interval.(r1[ind].coeffs, r2[ind].coeffs))
        @test all(isguaranteed.(getfield(r1[ind], :coeffs)))
        @test all(isguaranteed.(getfield(r2[ind], :coeffs)))
    end


    # Tests `evaluate`
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
        " ( Interval{Float64}(1.5, 2.0, com) + Interval{Float64}(0.0, 0.0, com)im ) +" *
        " ( Interval{Float64}(-2.0, -1.0, com) + Interval{Float64}(-1.0, 1.0, com)im ) t" *
        " + ( Interval{Float64}(-1.0, 1.5, com) + Interval{Float64}(-1.0, 1.5, com)im ) t²" *
        " + ( Interval{Float64}(0.0, 0.0, com) + Interval{Float64}(-1.0, 1.5, com)im ) t³ "
    displayBigO(true)
    @test string(Taylor1(vc, 5)) ==
        " ( Interval{Float64}(1.5, 2.0, com) + Interval{Float64}(0.0, 0.0, com)im ) +" *
        " ( Interval{Float64}(-2.0, -1.0, com) + Interval{Float64}(-1.0, 1.0, com)im ) t" *
        " + ( Interval{Float64}(-1.0, 1.5, com) + Interval{Float64}(-1.0, 1.5, com)im ) t²" *
        " + ( Interval{Float64}(0.0, 0.0, com) + Interval{Float64}(-1.0, 1.5, com)im ) t³ + 𝒪(t⁶)"


    # Iss 351 (inspired by a test in ReachabilityAnalysis)
    p1 = Taylor1([interval(0, 0), interval(0, 0.1) + interval(0, 0.01) * y], 4)
    p2 = Taylor1([interval(0, 0), interval(0, 0.5) + interval(0, 0.02) * x + interval(0, 0.03) * y], 4)
    @test evaluate([p1, p2], interval(0, 1)) == [p1[1], p2[1]]
    @test typeof(p1(interval(0, 1))) == TaylorN{Interval{Float64}}

    # Tests related to Iss #311
    # `sqrt` and `pow` defined on Interval(0,Inf)
    @test_throws DomainError sqrt(ti)
    @test sqrt(interval(0.0, 1.e-15) + ti) == sqrt(interval(-1.e-15, 1.e-15) + ti)
    aa = sqrt(sqrt(interval(0.0, 1.e-15) + ti))
    @test aa == sqrt(sqrt(interval(-1.e-15, 1.e-15) + ti))
    bb = (interval(0.0, 1.e-15) + ti)^(1/4)
    @test bb == (interval(-1.e-15, 1.e-15) + ti)^(1/4)
    @test all(isinterior.(aa.coeffs[2:end], bb.coeffs[2:end]))
    @test_throws DomainError sqrt(x)
    @test sqrt(interval(-1,1)+x) == sqrt(interval(0,1)+x)
    @test (interval(-1,1)+x)^(1/4) == (interval(0,1)+x)^(1/4)

    # `log` defined on Interval(0,Inf)
    @test_throws DomainError log(ti)
    @test log(interval(0.0, 1.e-15) + ti) == log(interval(-1.e-15, 1.e-15) + ti)
    @test_throws DomainError log(y)
    @test log(interval(0.0, 1.e-15) + y) == log(interval(-1.e-15, 1.e-15) + y)
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
