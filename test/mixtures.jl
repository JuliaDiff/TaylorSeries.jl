# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries
using Base.Test

@testset "Tests with mixures of Taylor1 and TaylorN" begin
    set_variables("x", numvars=2, order=6)
    xH = HomogeneousPolynomial(Int, 1)
    yH = HomogeneousPolynomial(Int, 2)
    tN = Taylor1(TaylorN{Float64}, 3)

    @test eltype(tN) == TaylorN{Float64}
    @test tN.order == 3
    @test string(zero(tN)) == "  0.0 + 𝒪(‖x‖¹) + 𝒪(t⁴)"
    @test string(tN) == " ( 1.0 + 𝒪(‖x‖¹)) t + 𝒪(t⁴)"
    @test string(tN + 3Taylor1(Int, 2)) == " ( 4.0 + 𝒪(‖x‖¹)) t + 𝒪(t⁴)"
    @test string(xH * tN) == " ( 1.0 x₁ + 𝒪(‖x‖²)) t + 𝒪(t⁴)"

    tN = Taylor1([zero(TaylorN(Float64,1)), one(TaylorN(Float64,1))], 3)
    @test typeof(tN) == Taylor1{TaylorN{Float64}}
    @test string(zero(tN)) == "  0.0 + 𝒪(‖x‖⁷) + 𝒪(t⁴)"
    @test string(tN) == " ( 1.0 + 𝒪(‖x‖⁷)) t + 𝒪(t⁴)"
    @test string(Taylor1([xH+yH])) == "  1 x₁ + 1 x₂ + 𝒪(t¹)"
    @test string(Taylor1([zero(xH), xH*yH])) == " ( 1 x₁ x₂) t + 𝒪(t²)"
    @test string(tN * Taylor1([0,TaylorN([xH+yH])])) ==
        " ( 1.0 x₁ + 1.0 x₂ + 𝒪(‖x‖⁷)) t² + 𝒪(t⁴)"

    t = Taylor1(3)
    xHt = HomogeneousPolynomial(typeof(t), 1)
    @test string(xHt) == " ( 1.0 + 𝒪(t¹)) x₁"
    xHt = HomogeneousPolynomial([one(t), zero(t)])
    yHt = HomogeneousPolynomial([zero(t), t])
    @test string(xHt) == " ( 1.0 + 𝒪(t⁴)) x₁"
    @test string(yHt) == " ( 1.0 t + 𝒪(t⁴)) x₂"
    @test string(HomogeneousPolynomial([t])) == " ( 1.0 t + 𝒪(t⁴))"
    @test 3*xHt == HomogeneousPolynomial([3*one(t), zero(t)])
    @test complex(0,1)*xHt == HomogeneousPolynomial([1im*one(t), zero(1im*t)])

    tN1 = TaylorN([HomogeneousPolynomial([t]),xHt,yHt^2])
    t1N = convert(Taylor1{TaylorN{Float64}}, tN1)
    ctN1 = convert(TaylorN{Taylor1{Float64}}, t1N)
    @test eltype(xHt) == Taylor1{Float64}
    @test eltype(tN1) == Taylor1{Float64}
    @test eltype(Taylor1([xH])) == HomogeneousPolynomial{Int64}
    @test eltype(tN1) == Taylor1{Float64}
    @test get_order(HomogeneousPolynomial([Taylor1(1), 1.0+Taylor1(2)])) == 1
    @test string(tN1) ==
        " ( 1.0 t + 𝒪(t⁴)) + ( 1.0 + 𝒪(t⁴)) x₁ + ( 1.0 t² + 𝒪(t⁴)) x₂² + 𝒪(‖x‖³)"
    @test string(t1N) ==
        "  1.0 x₁ + 𝒪(‖x‖³) + ( 1.0 + 𝒪(‖x‖³)) t + ( 1.0 x₂² + 𝒪(‖x‖³)) t² + 𝒪(t⁴)"
    @test tN1 == ctN1
    @test tN1+tN1 == 2*tN1
    @test tN1+1im*tN1 == complex(1,1)*tN1
    @test tN1+t == t+tN1
    @test tN1-t == -t+tN1
    @test tN1-tN1 == zero(tN1)
    @test string(t1N*t1N) ==
        "  1.0 x₁² + 𝒪(‖x‖³) + ( 2.0 x₁ + 𝒪(‖x‖³)) t + ( 1.0 + 𝒪(‖x‖³)) t² + ( 2.0 x₂² + 𝒪(‖x‖³)) t³ + 𝒪(t⁴)"
    @test !isnan(tN1)
    @test !isinf(tN1)

    @test mod(tN1+1,1.0) == 0+tN1
    @test mod(tN1-1.125,2) == 0.875+tN1
    @test (rem(tN1+1.125,1.0)).coeffs[1].coeffs[1] == 0.125 + t
    @test (rem(tN1-1.125,2)).coeffs[1].coeffs[1] == -1.125 + t
    @test mod2pi(-3pi+tN1).coeffs[1].coeffs[1].coeffs[1] ≈ pi
    @test mod2pi(0.125+2pi+tN1).coeffs[1].coeffs[1].coeffs[1] ≈ 0.125
    @test mod(t1N+1.125,1.0) == 0.125+t1N
    @test mod(t1N-1.125,2) == 0.875+t1N
    @test (rem(t1N+1.125,1.0)).coeffs[1] == 0.125 + t1N.coeffs[1]
    @test (rem(t1N-1.125,2)).coeffs[1] == -1.125 + t1N.coeffs[1]
    @test mod2pi(-3pi+t1N).coeffs[1].coeffs[1].coeffs[1] ≈ pi
    @test mod2pi(0.125+2pi+t1N).coeffs[1].coeffs[1].coeffs[1] ≈ 0.125

    @test convert(Array{Taylor1{TaylorN{Float64}},1}, [tN1, tN1]) == [t1N, t1N]
    @test convert(Array{Taylor1{TaylorN{Float64}},2}, [tN1 tN1]) == [t1N t1N]
    @test convert(Array{TaylorN{Taylor1{Float64}},1}, [t1N, t1N]) == [tN1, tN1]
    @test convert(Array{TaylorN{Taylor1{Float64}},2}, [t1N t1N]) == [tN1 tN1]

    @test string(evaluate(t1N, 0.0)) == " 1.0 x₁ + 𝒪(‖x‖³)"
    @test string(evaluate(t1N^2, 1.0)) == " 1.0 + 2.0 x₁ + 1.0 x₁² + 2.0 x₂² + 𝒪(‖x‖³)"
    v = zeros(TaylorN{Float64},2)
    @test evaluate!([t1N, t1N^2], 0.0, v) == nothing
    @test v[1] == TaylorN([xHt])
    @test v[2] == TaylorN([xHt^2])

    tint = Taylor1(Int, 10)
    t = Taylor1(10)
    x = TaylorN( [HomogeneousPolynomial(zero(t)), HomogeneousPolynomial([one(t),zero(t)])], 5)
    y = TaylorN(typeof(tint), 2, order=5)
    @test typeof(x) == TaylorN{Taylor1{Float64}}
    @test eltype(y) == Taylor1{Int}
    @test -x == 0 - x
    @test +y == y
    @test one(y)/(1+x) == 1 - x + x^2 - x^3 + x^4 - x^5
    @test one(y)/(1+y) == 1 - y + y^2 - y^3 + y^4 - y^5

    # See #92 and #94
    δx, δy = set_variables("δx δy")
    xx = 1+Taylor1(δx,5)
    @test typeof(xx) == Taylor1{TaylorN{Float64}}
    @test eltype(xx) == TaylorN{Float64}
    @test !isnan(xx)
    @test !isnan(δx)
    @test !isinf(xx)
    @test !isinf(δx)
    @test +xx == xx
    @test -xx == 0 - xx
    @test xx/1.0 == 1.0*xx
    @test xx + xx == xx*2
    @test xx - xx == zero(xx)
    @test xx*xx == xx^2
    @test xx/xx == one(xx)
    @test xx/(1+δx) == one(xx)
end
