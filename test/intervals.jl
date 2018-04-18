# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries, IntervalArithmetic

if VERSION < v"0.7.0-DEV.2004"
    using Base.Test
    eeuler = Base.e
else
    using Test
    eeuler = Base.MathConstants.e
end

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
    @test all((p5(x,-b))[:] .⊆ ((x-b)^5)[:])

end
