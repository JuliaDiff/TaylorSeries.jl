# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries
using Base.Test

@testset "Tests with mixtures of Taylor1 and TaylorN" begin
    @test TaylorSeries.NumberNotSeries == Union{Real,Complex}
    @test TaylorSeries.NumberNotSeriesN == Union{Real,Complex,Taylor1}

    set_variables("x", numvars=2, order=6)
    xH = HomogeneousPolynomial(Int, 1)
    yH = HomogeneousPolynomial(Int, 2)
    tN = Taylor1(TaylorN{Float64}, 3)

    @test eltype(xH) == Int
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
    @test eltype(xHt) == Taylor1{Float64}
    @test string(xHt) == " ( 1.0 + 𝒪(t¹)) x₁"
    xHt = HomogeneousPolynomial([one(t), zero(t)])
    yHt = HomogeneousPolynomial([zero(t), t])
    @test string(xHt) == " ( 1.0 + 𝒪(t⁴)) x₁"
    @test string(yHt) == " ( 1.0 t + 𝒪(t⁴)) x₂"
    @test string(HomogeneousPolynomial([t])) == " ( 1.0 t + 𝒪(t⁴))"
    @test 3*xHt == HomogeneousPolynomial([3*one(t), zero(t)])
    @test t*xHt == HomogeneousPolynomial([t, zero(t)])
    @test complex(0,1)*xHt == HomogeneousPolynomial([1im*one(t), zero(1im*t)])
    @test eltype(complex(0,1)*xHt) == Taylor1{Complex128}

    tN1 = TaylorN([HomogeneousPolynomial([t]),xHt,yHt^2])
    @test tN1[1] == HomogeneousPolynomial([t])
    t1N = convert(Taylor1{TaylorN{Float64}}, tN1)
    @test t1N[1] == HomogeneousPolynomial(1)
    ctN1 = convert(TaylorN{Taylor1{Float64}}, t1N)
    @test eltype(xHt) == Taylor1{Float64}
    @test eltype(tN1) == Taylor1{Float64}
    @test eltype(Taylor1([xH])) == HomogeneousPolynomial{Int64}
    @test eltype(tN1) == Taylor1{Float64}
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
    @test tN1-tN1 == zero(tN1)
    @test string(t1N*t1N) ==
        "  1.0 x₁² + 𝒪(‖x‖³) + ( 2.0 x₁ + 𝒪(‖x‖³)) t + ( 1.0 + 𝒪(‖x‖³)) t² + ( 2.0 x₂² + 𝒪(‖x‖³)) t³ + 𝒪(t⁴)"
    @test !isnan(tN1)
    @test !isinf(tN1)

    @test mod(tN1+1,1.0) == 0+tN1
    @test mod(tN1-1.125,2) == 0.875+tN1
    @test (rem(tN1+1.125,1.0))[1][1] == 0.125 + t
    @test (rem(tN1-1.125,2))[1][1] == -1.125 + t
    @test mod2pi(-3pi+tN1)[1][1][1] ≈ pi
    @test mod2pi(0.125+2pi+tN1)[1][1][1] ≈ 0.125
    @test mod(t1N+1.125,1.0) == 0.125+t1N
    @test mod(t1N-1.125,2) == 0.875+t1N
    @test (rem(t1N+1.125,1.0))[1] == 0.125 + t1N[1]
    @test (rem(t1N-1.125,2))[1] == -1.125 + t1N[1]
    @test mod2pi(-3pi+t1N)[1][1][1] ≈ pi
    @test mod2pi(0.125+2pi+t1N)[1][1][1] ≈ 0.125

    @test abs(tN1+1) == 1+tN1
    @test abs(tN1-1) == 1-tN1
    @test_throws ArgumentError abs(tN1)
    @test_throws ArgumentError abs(t1N)

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
    @test (1+y)/one(t) == 1 + y
    @test typeof(y+t) == TaylorN{Taylor1{Float64}}

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
    @test xx*δx + Taylor1(typeof(δx),5) == δx + δx^2 + Taylor1(typeof(δx),5)
    @test xx/(1+δx) == one(xx)
    @test typeof(xx+δx) == Taylor1{TaylorN{Float64}}

    X, Y = set_variables(Taylor1{Float64}, "x y")
    @test typeof( norm(X) ) == Float64
    @test norm(X) > 0
    @test norm(X+Y) == sqrt(2)
    @test norm(-10X+4Y,Inf) == 10.

    @test TaylorSeries.rtoldefault(TaylorN{Taylor1{Int64}}) == 0
    @test TaylorSeries.rtoldefault(Taylor1{TaylorN{Int64}}) == 0
    for T in (Float64, BigFloat)
        @test TaylorSeries.rtoldefault(TaylorN{Taylor1{T}}) == sqrt(eps(T))
        @test TaylorSeries.rtoldefault(Taylor1{TaylorN{T}}) == sqrt(eps(T))
        @test TaylorSeries.real(TaylorN{Taylor1{T}}) == TaylorN{Taylor1{T}}
        @test TaylorSeries.real(Taylor1{TaylorN{T}}) == Taylor1{TaylorN{T}}
        @test TaylorSeries.real(TaylorN{Taylor1{Complex{T}}}) == TaylorN{Taylor1{T}}
        @test TaylorSeries.real(Taylor1{TaylorN{Complex{T}}}) == Taylor1{TaylorN{T}}
    end

    rndT1(ord1) = Taylor1(-1+2rand(ord1+1)) # generates a random Taylor1 with order `ord`
    nmonod(s, d) = binomial(d+s-1, d) #number of monomials in s variables with exact degree d
    #rndHP generates a random `ordHP`-th order homog. pol. of Taylor1s, each with order `ord1`
    rndHP(ordHP, ord1) = HomogeneousPolynomial( [rndT1(ord1) for i in 1:nmonod(get_numvars(), ordHP)] )
    #rndTN generates a random `ordHP`-th order TaylorN of of Taylor1s, each with order `ord1`
    rndTN(ordN, ord1) = TaylorN([rndHP(i, ord1) for i in 0:ordN])

    P = rndTN(get_order(), 3)
    @test P ≈ P
    Q = deepcopy(P)
    Q[2][2] = Taylor1([NaN, Inf])
    @test isnan(Q)
    @test isinf(Q)
    @test !isfinite(Q)
    Q[2][2] = P[2][2]+sqrt(eps())/2
    @test isapprox(P, Q, rtol=1.0)
    Q[2][2] = P[2][2]+10sqrt(eps())
    @test !isapprox(P, Q, atol=sqrt(eps()), rtol=0)
    Q[2][2] = P[2][2]+eps()/2
    @test isapprox(Q, Q, atol=eps(), rtol=0)
    @test isapprox(Q, P, atol=eps(), rtol=0)
    Q[2][1] = P[2][1]-10eps()
    @test !isapprox(Q, P, atol=eps(), rtol=0)

    p1N = Taylor1([X^2,X*Y,Y+X,Y^2])
    q1N = Taylor1([X^2,(1.0+sqrt(eps()))*X*Y,Y+X,Y^2])
    @test p1N ≈ p1N
    @test p1N ≈ q1N
end
