# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries

using Test
using LinearAlgebra, SparseArrays

# This is used to check the fallack of pretty_print
struct SymbNumber <: Number
    s :: Symbol
end
Base.iszero(::SymbNumber) = false

@testset "Tests for Taylor1 expansions" begin
    eeuler = Base.MathConstants.e

    ta(a) = Taylor1([a,one(a)],15)
    t = Taylor1(Int,15)
    tim = im*t
    zt = zero(t)
    ot = 1.0*one(t)
    tol1 = eps(1.0)

    @test TS === TaylorSeries
    @test Taylor1 <: AbstractSeries
    @test Taylor1{Float64} <: AbstractSeries{Float64}
    @test TS.numtype(1.0) == eltype(1.0)
    @test TS.numtype([1.0]) == eltype([1.0])
    @test TS.normalize_taylor(t) == t
    @test TS.normalize_taylor(tim) == tim

    @test Taylor1([1,2,3,4,5], 2) == Taylor1([1,2,3])
    @test Taylor1(t[0:3]) == Taylor1(t[0:get_order(t)], 4)
    @test get_order(t) == 15
    @test get_order(Taylor1([1,2,3,4,5], 2)) == 2

    @test size(t) == (16,)
    @test firstindex(t) == 0
    @test lastindex(t) == 15
    @test eachindex(t) == 0:15
    @test iterate(t) == (0.0, 1)
    @test iterate(t, 1) == (1.0, 2)
    @test iterate(t, 16) == nothing
    @test axes(t) == ()
    @test axes([t]) == (Base.OneTo(1),)

    @testset "Total order" begin
        @test 1 + t ≥ 1.0 > 0.5 + t > t^2 ≥ zero(t)
        @test -1.0 < -1/1000 - t  < -t < -t^2 ≤ 0
    end

    v = [1,2]
    @test typeof(TS.resize_coeffs1!(v,3)) == Nothing
    @test v == [1,2,0,0]
    TS.resize_coeffs1!(v,0)
    @test v == [1]
    TS.resize_coeffs1!(v,3)
    setindex!(Taylor1(v),3,2)
    @test v == [1,0,3,0]
    pol_int = Taylor1(v)
    @test pol_int[:] == [1,0,3,0]
    @test pol_int[:] == pol_int.coeffs[:]
    @test pol_int[1:2:3] == pol_int.coeffs[2:2:4]
    setindex!(pol_int,0,0:2)
    @test v == zero(v)
    setindex!(pol_int,1,:)
    @test v == ones(Int, 4)
    setindex!(pol_int, v, :)
    @test v == ones(Int, 4)
    setindex!(pol_int, zeros(Int, 4), 0:3)
    @test v == zeros(Int, 4)
    pol_int[:] .= 0
    @test v == zero(v)
    pol_int[0:2:end] = 2
    @test all(v[1:2:end] .== 2)
    pol_int[0:2:3] = [0, 1]
    @test all(v[1:2:3] .== [0, 1])
    rv = [rand(0:3) for i in 1:4]
    @test Taylor1(rv)[:] == rv
    y = sin(Taylor1(16))
    @test y[:] == y.coeffs
    y[:] .= cos(Taylor1(16))[:]
    @test y == cos(Taylor1(16))
    @test y[:] == cos(Taylor1(16))[:]
    y = sin(Taylor1(16))
    rv = rand(5)
    y[0:4] .= rv
    @test y[0:4] == rv
    @test y[5:end] == y.coeffs[6:end]
    rv = rand( length(y.coeffs) )
    y[:] .= rv
    @test y[:] == rv
    y[:] .= cos(Taylor1(16)).coeffs
    @test y == cos(Taylor1(16))
    @test y[:] == cos(Taylor1(16))[:]
    y[:] .= 0.0
    @test y[:] == zero(y[:])
    y = sin(Taylor1(16))
    rv = rand.(length(0:4))
    y[0:4] .= rv
    @test y[0:4] == rv
    @test y[6:end] == sin(Taylor1(16))[6:end]
    rv = rand.(length(y))
    y[:] .= rv
    @test y[:] == rv
    y[0:4:end] .= 1.0
    @test all(y.coeffs[1:4:end] .== 1.0)
    y[0:2:8] .= rv[1:5]
    @test y.coeffs[1:2:9] == rv[1:5]
    @test_throws AssertionError y[0:2:3] = rv

    @test Taylor1([0,1,0,0]) == Taylor1(3)
    @test getcoeff(Taylor1(Complex{Float64},3),1) == complex(1.0,0.0)
    @test Taylor1(Complex{Float64},3)[1] == complex(1.0,0.0)
    @test getindex(Taylor1(3),1) == 1.0
    @inferred convert(Taylor1{Complex{Float64}},ot) == Taylor1{Complex{Float64}}
    @test eltype(convert(Taylor1{Complex{Float64}},ot)) == Taylor1{Complex{Float64}}
    @test eltype(convert(Taylor1{Complex{Float64}},1)) == Taylor1{Complex{Float64}}
    @test eltype(convert(Taylor1, 1im)) == Taylor1{Complex{Int}}
    @test TS.numtype(convert(Taylor1{Complex{Float64}},ot)) == Complex{Float64}
    @test TS.numtype(convert(Taylor1{Complex{Float64}},1)) == Complex{Float64}
    @test TS.numtype(convert(Taylor1, 1im)) == Complex{Int}
    @test convert(Taylor1, 1im) == Taylor1(1im, 0)
    @test convert(eltype(t), t) === t
    @test convert(eltype(ot), ot) === ot
    @test convert(Taylor1{Int},[0,2]) == 2*t
    @test convert(Taylor1{Complex{Int}},[0,2]) == (2+0im)*t
    @test convert(Taylor1{BigFloat},[0.0, 1.0]) == ta(big(0.0))
    @test promote(t,Taylor1(1.0,0)) == (ta(0.0),ot)
    @test promote(0,Taylor1(1.0,0)) == (zt,ot)
    @test eltype(promote(ta(0),zeros(Int,2))[2]) == Taylor1{Int}
    @test eltype(promote(ta(0.0),zeros(Int,2))[2]) == Taylor1{Float64}
    @test eltype(promote(0,Taylor1(ot))[1]) == Taylor1{Float64}
    @test eltype(promote(1.0+im, zt)[1]) == Taylor1{Complex{Float64}}
    @test TS.numtype(promote(ta(0),zeros(Int,2))[2]) == Int
    @test TS.numtype(promote(ta(0.0),zeros(Int,2))[2]) == Float64
    @test TS.numtype(promote(0,Taylor1(ot))[1]) == Float64
    @test TS.numtype(promote(1.0+im, zt)[1]) == Complex{Float64}

    @test length(Taylor1(10)) == 11
    @test length.( TS.fixorder(zt, Taylor1([1])) ) == (16, 16)
    @test length.( TS.fixorder(zt, Taylor1([1], 1)) ) == (2, 2)
    @test eltype(TS.fixorder(zt,Taylor1([1]))[1]) == Taylor1{Int}
    @test TS.numtype(TS.fixorder(zt,Taylor1([1]))[1]) == Int
    @test findfirst(t) == 1
    @test findfirst(t^2) == 2
    @test findfirst(ot) == 0
    @test findfirst(zt) == -1
    @test iszero(zero(t))
    @test !iszero(one(t))
    @test @inferred isinf(Taylor1([typemax(1.0)]))
    @test @inferred isnan(Taylor1([typemax(1.0), NaN]))

    @test constant_term(2.0) == 2.0
    @test constant_term(t) == 0
    @test constant_term(tim) == complex(0, 0)
    @test constant_term([zt, t]) == [0, 0]
    @test linear_polynomial(2) == 2
    @test linear_polynomial(t) == t
    @test linear_polynomial(1+tim^2) == zero(tim)
    @test get_order(linear_polynomial(1+tim^2)) == get_order(tim)
    @test linear_polynomial([zero(tim), tim, tim^2]) == [zero(tim), tim, zero(tim)]
    @test nonlinear_polynomial(2im) == 0im
    @test nonlinear_polynomial(1+t) == zero(t)
    @test nonlinear_polynomial(1+tim^2) == tim^2
    @test nonlinear_polynomial([zero(tim), tim, 1+tim^2]) == [zero(tim), zero(tim), tim^2]

    @test ot == 1
    @test 0.0 == zt
    @test getcoeff(tim,1) == complex(0,1)
    @test zt+1.0 == ot
    @test 1.0-ot == zt
    @test t+t == 2t
    @test t-t == zt
    @test +t == -(-t)

    tsquare = Taylor1([0,0,1],15)
    @test t * true == t
    @test false * t == zero(t)
    @test t^0 == t^0.0 == one(t)
    @test t*t == tsquare
    @test t*1 == t
    @test 0*t == zt
    @test (-t)^2 == tsquare
    @test t^3 == tsquare*t
    @test zero(t)/t == zero(t)
    @test get_order(zero(t)/t) == get_order(t)
    @test one(t)/one(t) == 1.0
    @test tsquare/t == t
    @test get_order(tsquare/t) == get_order(tsquare)-1
    @test t/(t*3) == (1/3)*ot
    @test get_order(t/(t*3)) == get_order(t)-1
    @test t/3im == -tim/3
    @test 1/(1-t) == Taylor1(ones(t.order+1))
    @test Taylor1([0,1,1])/t == t+1
    @test get_order(Taylor1([0,1,1])/t) == 1
    @test (t+im)^2 == tsquare+2im*t-1
    @test (t+im)^3 == Taylor1([-1im,-3,3im,1],15)
    @test (t+im)^4 == Taylor1([1,-4im,-6,4im,1],15)
    @test imag(tsquare+2im*t-1) == 2t
    @test (Rational(1,2)*tsquare)[2] == 1//2
    @test t^2/tsquare == ot
    @test get_order(t^2/tsquare) == get_order(t)-2
    @test ((1+t)^(1/3))[2]+1/9 ≤ tol1
    @test (1.0-tsquare)^3 == (1.0-t)^3*(1.0+t)^3
    @test (1-tsquare)^2 == (1+t)^2.0 * (1-t)^2.0
    @test (sqrt(1+t))[2] == -1/8
    @test ((1-tsquare)^(1//2))^2 == 1-tsquare
    @test ((1-t)^(1//4))[14] == -4188908511//549755813888
    @test abs(((1+t)^3.2)[13] + 5.4021062656e-5) < tol1
    @test Taylor1(BigFloat,5)/6 == 1im*Taylor1(5)/complex(0,BigInt(6))
    @test Taylor1(BigFloat,5)/(6*Taylor1(3)) == 1/BigInt(6)
    @test Taylor1(BigFloat,5)/(6im*Taylor1(3)) == -1im/BigInt(6)
    @test isapprox((1+(1.5+t)/4)^(-2), inv(1+(1.5+t)/4)^2, rtol=eps(Float64))
    @test isapprox((1+(big(1.5)+t)/4)^(-2), inv(1+(big(1.5)+t)/4)^2, rtol=eps(BigFloat))
    @test isapprox((1+(1.5+t)/4)^(-2), inv(1+(1.5+t)/4)^2, rtol=eps(Float64))
    @test isapprox((1+(big(1.5)+t)/4)^(-2), inv(1+(big(1.5)+t)/4)^2, rtol=eps(BigFloat))
    @test isapprox((1+(1.5+t)/5)^(-2.5), inv(1+(1.5+t)/5)^2.5, rtol=eps(Float64))
    @test isapprox((1+(big(1.5)+t)/5)^(-2.5), inv(1+(big(1.5)+t)/5)^2.5, rtol=2eps(BigFloat))

    # These tests involve some sort of factorization
    @test t/(t+t^2) == 1/(1+t)
    @test get_order(t/(t+t^2)) == get_order(1/(1+t))-1
    @test sqrt(t^2+t^3) == t*sqrt(1+t)
    @test get_order(sqrt(t^2+t^3)) == get_order(t) >> 1
    @test get_order(t*sqrt(1+t)) == get_order(t)
    @test (t^3+t^4)^(1/3) ≈ t*(1+t)^(1/3)
    @test norm((t^3+t^4)^(1/3) - t*(1+t)^(1/3), Inf) < eps()
    @test get_order((t^3+t^4)^(1/3)) == 5
    @test ((t^3+t^4)^(1/3))[5] == -10/243

    trational = ta(0//1)
    @inferred ta(0//1) == Taylor1{Rational{Int}}
    @test eltype(trational) == Taylor1{Rational{Int}}
    @test TS.numtype(trational) == Rational{Int}
    @test trational + 1//3 == Taylor1([1//3,1],15)
    @test complex(3,1)*trational^2 == Taylor1([0//1,0//1,complex(3,1)//1],15)
    @test trational^2/3 == Taylor1([0//1,0//1,1//3],15)
    @test trational^3/complex(7,1) == Taylor1([0,0,0,complex(7//50,-1//50)],15)
    @test sqrt(zero(t)) == zero(t)

    @test isapprox( rem(4.1 + t,4)[0], 0.1 )
    @test isapprox( mod(4.1 + t,4)[0], 0.1 )
    @test isapprox( rem(1+Taylor1(Int,4),4.0)[0], 1.0 )
    @test isapprox( mod(1+Taylor1(Int,4),4.0)[0], 1.0 )
    @test isapprox( mod2pi(2pi+0.1+t)[0], 0.1 )

    @test abs(ta(1)) == ta(1)
    @test abs(ta(-1.0)) == -ta(-1.0)


    @test taylor_expand(x->2x,order=10) == 2*Taylor1(10)
    @test taylor_expand(x->x^2+1) == Taylor1(15)*Taylor1(15) + 1
    @test evaluate(taylor_expand(cos,0.)) == cos(0.)
    @test evaluate(taylor_expand(tan,pi/4)) == tan(pi/4)
    @test eltype(taylor_expand(x->x^2+1,1)) == Taylor1{Int}
    @test TS.numtype(taylor_expand(x->x^2+1,1)) == Int
    tsq = t^2
    update!(tsq,2.0)
    @test tsq == (t+2.0)^2
    update!(tsq,-2)
    @test tsq == t^2

    @test log(exp(tsquare)) == tsquare
    @test exp(log(1-tsquare)) == 1-tsquare
    @test constant_term(expm1(1.0e-16+t)) == 1.0e-16
    @test expm1(1.e-16+t).coeffs[2:end] ≈ expm1(1.e-16).*(exp(t)-1).coeffs[2:end]
    @test log((1-t)^2) == 2*log(1-t)
    @test log1p(0.25 + t) == log(1.25+t)
    @test log1p(-t^2) == log(1-t^2)

    st, ct = sincos(t)
    @test real(exp(tim)) == ct
    @test imag(exp(tim)) == st
    @test exp(conj(tim)) == ct-im*st == exp(tim')
    st, ct = sincospi(t)
    @test (st, ct) == sincos(pi*t)
    @test real(exp(pi*tim)) == cospi(t)
    @test imag(exp(pi*tim)) == sinpi(t)
    @test exp(pi*conj(tim)) == ct-im*st == exp(pi*tim')
    @test abs2(tim) == tsquare
    @test abs(tim) == t
    @test isapprox(abs2(exp(tim)), ot)
    @test isapprox(abs(exp(tim)), ot)
    @test (exp(t))^(2im) == cos(2t)+im*sin(2t)
    @test (exp(t))^Taylor1([-5.2im]) == cos(5.2t)-im*sin(5.2t)
    @test getcoeff(convert(Taylor1{Rational{Int}},cos(t)),8) == 1//factorial(8)
    @test abs((tan(t))[7]- 17/315) < tol1
    @test abs((tan(t))[13]- 21844/6081075) < tol1
    @test tan(1.3+t) ≈ sin(1.3+t)/cos(1.3+t)
    @test cot(1.3+t) ≈ 1/tan(1.3+t)
    @test evaluate(exp(Taylor1([0,1],17)),1.0) == 1.0*eeuler
    @test evaluate(exp(Taylor1([0,1],1))) == 1.0
    @test evaluate(exp(t),t^2) == exp(t^2)
    @test evaluate(exp(Taylor1(BigFloat, 15)), t^2) == exp(Taylor1(BigFloat, 15)^2)
    @test evaluate(exp(Taylor1(BigFloat, 15)), t^2) isa Taylor1{BigFloat}
    #Test function-like behavior for Taylor1s
    t17 = Taylor1([0,1],17)
    myexpfun = exp(t17)
    @test myexpfun(1.0) == 1.0*eeuler
    @test myexpfun() == 1.0
    @test myexpfun(t17^2) == exp(t17^2)
    @test exp(t17^2)(t17) == exp(t17^2)
    q, p = sincospi(t17)
    @test cospi(-im*t)(1)+im*sinpi(-im*t)(1) == exp(-im*pi*t)(im)
    @test p(-im*t17)(1)+im*q(-im*t17)(1) ≈ exp(-im*pi*t17)(im)
    q, p = sincos(t17)
    @test cos(-im*t)(1)+im*sin(-im*t)(1) == exp(-im*t)(im)
    @test p(-im*t17)(1)+im*q(-im*t17)(1) == exp(-im*t17)(im)
    cossin1 = x->p(q(x))
    @test evaluate(p, evaluate(q, pi/4)) == cossin1(pi/4)
    cossin2 = p(q)
    @test evaluate(evaluate(p,q), pi/4) == cossin2(pi/4)
    @test evaluate(p, q) == cossin2
    @test p(q)() == evaluate(evaluate(p, q))
    @test evaluate(p, q) == p(q)
    @test evaluate(q, p) == q(p)
    cs = x->cos(sin(x))
    csdiff = (cs(t17)-cossin2(t17)).(-2:0.1:2)
    @test norm(csdiff, 1) < 5e-15
    a = [p,q]
    @test a(0.1) == evaluate.([p,q],0.1)
    @test a.(0.1) == a(0.1)
    @test a.() == evaluate.([p, q])
    @test a.() == [p(), q()]
    @test a.() == a()
    @test view(a, 1:1)() == [a[1]()]
    vr = rand(2)
    @test p.(vr) == evaluate.([p], vr)
    Mr = rand(3,3,3)
    @test p.(Mr) == evaluate.([p], Mr)
    vr = rand(5)
    @test p(vr) == p.(vr)
    @test view(a, 1:1)(vr) == evaluate.([p],vr)
    @test p(Mr) == p.(Mr)
    @test p(Mr) == evaluate.([p], Mr)
    taylor_a = Taylor1(Int,10)
    taylor_x = exp(Taylor1(Float64,13))
    @test taylor_x(taylor_a) == evaluate(taylor_x, taylor_a)
    A_T1 = [t 2t 3t; 4t 5t 6t ]
    @test evaluate(A_T1,1.0) == [1.0  2.0  3.0; 4.0  5.0  6.0]
    @test evaluate(A_T1,1.0) == A_T1(1.0)
    @test evaluate(A_T1) == A_T1()
    @test A_T1(tsquare) == [tsquare 2tsquare 3tsquare; 4tsquare 5tsquare 6tsquare]
    @test view(A_T1, :, :)(1.0) == A_T1(1.0)
    @test view(A_T1, :, 1)(1.0) == A_T1[:,1](1.0)

    @test sin(asin(tsquare)) == tsquare
    @test tan(atan(tsquare)) == tsquare
    @test atan(tan(tsquare)) == tsquare
    @test atan(sin(tsquare)/cos(tsquare)) == atan(sin(tsquare), cos(tsquare))
    @test constant_term(atan(sin(3pi/4+tsquare), cos(3pi/4+tsquare))) == 3pi/4
    @test atan(sin(3pi/4+tsquare)/cos(3pi/4+tsquare)) - atan(sin(3pi/4+tsquare), cos(3pi/4+tsquare)) == -pi

    @test sinh(asinh(tsquare)) ≈ tsquare
    @test tanh(atanh(tsquare)) ≈ tsquare
    @test atanh(tanh(tsquare)) ≈ tsquare
    @test asinh(t) ≈ log(t + sqrt(t^2 + 1))
    @test cosh(asinh(t)) ≈ sqrt(t^2 + 1)

    t_complex = Taylor1(Complex{Int}, 15) # for use with acosh, which in the Reals is only defined for x ≥ 1
    @test cosh(acosh(t_complex)) ≈ t_complex
    @test differentiate(acosh(t_complex)) ≈ 1/sqrt(t_complex^2 - 1)
    @test acosh(t_complex) ≈ log(t_complex + sqrt(t_complex^2 - 1))
    @test sinh(acosh(t_complex)) ≈ sqrt(t_complex^2 - 1)

    @test asin(t) + acos(t) == pi/2
    @test differentiate(acos(t)) == - 1/sqrt(1-Taylor1(t.order-1)^2)
    @test get_order(differentiate(acos(t))) == t.order-1

    @test - sinh(t) + cosh(t) == exp(-t)
    @test  sinh(t) + cosh(t) == exp(t)
    @test evaluate(- sinh(t)^2 + cosh(t)^2 , rand()) == 1
    @test evaluate(- sinh(t)^2 + cosh(t)^2 , 0) == 1
    @test tanh(t + 0im) == -1im * tan(t*1im)
    @test  evaluate(tanh(t/2),1.5) == evaluate(sinh(t) / (cosh(t) + 1),1.5)
    @test cosh(t) == real(cos(im*t))
    @test sinh(t) == imag(sin(im*t))

    ut = 1.0*t
    tt = zero(ut)
    TS.one!(tt, ut, 0)
    @test tt[0] == 1.0
    TS.one!(tt, ut, 1)
    @test tt[1] == 0.0
    TS.abs!(tt, 1.0+ut, 0)
    @test tt[0] == 1.0
    TS.add!(tt, ut, ut, 1)
    @test tt[1] == 2.0
    TS.add!(tt, -3.0, 0)
    @test tt[0] == -3.0
    TS.add!(tt, -3.0, 1)
    @test tt[1] == 0.0
    TS.subst!(tt, ut, ut, 1)
    @test tt[1] == 0.0
    TS.subst!(tt, -3.0, 0)
    @test tt[0] == 3.0
    TS.subst!(tt, -2.5, 1)
    @test tt[1] == 0.0
    iind, cind = TS.divfactorization(ut, ut)
    @test iind == 1
    @test cind == 1.0
    TS.div!(tt, ut, ut, 0)
    @test tt[0] == cind
    TS.div!(tt, 1+ut, 1+ut, 0)
    @test tt[0] == 1.0
    TS.div!(tt, 1, 1+ut, 0)
    @test tt[0] == 1.0
    aux = tt
    TS.pow!(tt, 1.0+t, aux, 1.5, 0)
    @test tt[0] == 1.0
    TS.pow!(tt, 0.0*t, aux, 1.5, 0)
    @test tt[0] == 0.0
    TS.pow!(tt, 0.0+t, aux, 18, 0)
    @test tt[0] == 0.0
    TS.pow!(tt, 1.0+t, aux, 1.5, 0)
    @test tt[0] == 1.0
    TS.pow!(tt, 1.0+t, aux, 0.5, 1)
    @test tt[1] == 0.5
    TS.pow!(tt, 1.0+t, aux, 0, 0)
    @test tt[0] == 1.0
    TS.pow!(tt, 1.0+t, aux, 1, 1)
    @test tt[1] == 1.0
    tt = zero(ut)
    aux = tt
    TS.pow!(tt, 1.0+t, aux, 2, 0)
    @test tt[0] == 1.0
    TS.pow!(tt, 1.0+t, aux, 2, 1)
    @test tt[1] == 2.0
    TS.pow!(tt, 1.0+t, aux, 2, 2)
    @test tt[2] == 1.0
    TS.sqrt!(tt, 1.0+t, 0, 0)
    @test tt[0] == 1.0
    TS.sqrt!(tt, 1.0+t, 0)
    @test tt[0] == 1.0
    TS.exp!(tt, 1.0*t, 0)
    @test tt[0] == exp(t[0])
    TS.log!(tt, 1.0+t, 0)
    @test tt[0] == 0.0
    TS.log1p!(tt, 0.25+t, 0)
    @test tt[0] == log1p(0.25)
    TS.log1p!(tt, 0.25+t, 1)
    @test tt[1] == 1/1.25
    ct = zero(ut)
    TS.sincos!(tt, ct, 1.0*t, 0)
    @test tt[0] == sin(t[0])
    @test ct[0] == cos(t[0])
    TS.tan!(tt, 1.0*t, ct, 0)
    @test tt[0] == tan(t[0])
    @test ct[0] == tan(t[0])^2
    TS.asin!(tt, 1.0*t, ct, 0)
    @test tt[0] == asin(t[0])
    @test ct[0] == sqrt(1.0-t[0]^2)
    TS.acos!(tt, 1.0*t, ct, 0)
    @test tt[0] == acos(t[0])
    @test ct[0] == sqrt(1.0-t[0]^2)
    TS.atan!(tt, ut, ct, 0)
    @test tt[0] == atan(t[0])
    @test ct[0] == 1.0+t[0]^2
    TS.sinhcosh!(tt, ct, ut, 0)
    @test tt[0] == sinh(t[0])
    @test ct[0] == cosh(t[0])
    TS.tanh!(tt, ut, ct, 0)
    @test tt[0] == tanh(t[0])
    @test ct[0] == tanh(t[0])^2

    v = [sin(t), exp(-t)]
    vv = Vector{Float64}(undef, 2)
    @test evaluate!(v, zero(Int), vv) == nothing
    @test vv == [0.0,1.0]
    @test evaluate!(v, 0.0, vv) == nothing
    @test vv == [0.0,1.0]
    @test evaluate!(v, 0.0, view(vv, 1:2)) == nothing
    @test vv == [0.0,1.0]
    @test evaluate(v) == vv
    @test isapprox(evaluate(v, complex(0.0,0.2)),
        [complex(0.0,sinh(0.2)),complex(cos(0.2),sin(-0.2))], atol=eps(), rtol=0.0)
    m = [sin(t) exp(-t); cos(t) exp(t)]
    m0 = 0.5
    mres = Matrix{Float64}(undef, 2, 2)
    mres_expected = [sin(m0) exp(-m0); cos(m0) exp(m0)]
    @test evaluate!(m, m0, mres) == nothing
    @test mres == mres_expected

    ee_ta = exp(ta(1.0))
    @test get_order(differentiate(ee_ta, 0)) == 15
    @test get_order(differentiate(ee_ta, 1)) == 14
    @test get_order(differentiate(ee_ta, 16)) == 0
    @test differentiate(ee_ta, 0) == ee_ta
    expected_result_approx = Taylor1(ee_ta[0:10])
    @test differentiate(exp(ta(1.0)), 5) ≈ expected_result_approx atol=eps() rtol=0.0
    expected_result_approx = Taylor1(zero(ee_ta),0)
    @test differentiate(ee_ta, 16) == Taylor1(zero(ee_ta),0)
    @test eltype(differentiate(ee_ta, 16)) == eltype(ee_ta)
    ee_ta = exp(ta(1.0pi))
    expected_result_approx = Taylor1(ee_ta[0:12])
    @test differentiate(ee_ta, 3) ≈ expected_result_approx atol=eps(16.0) rtol=0.0
    expected_result_approx = Taylor1(ee_ta[0:5])
    @test differentiate(exp(ta(1.0pi)), 10) ≈ expected_result_approx atol=eps(64.0) rtol=0.0



    @test differentiate(exp(ta(1.0)), 5)() == exp(1.0)
    @test differentiate(exp(ta(1.0pi)), 3)() == exp(1.0pi)
    @test isapprox(derivative(exp(ta(1.0pi)), 10)() , exp(1.0pi) )

    @test differentiate(5, exp(ta(1.0))) == exp(1.0)
    @test differentiate(3, exp(ta(1.0pi))) == exp(1.0pi)
    @test isapprox(differentiate(10, exp(ta(1.0pi))) , exp(1.0pi) )

    @test integrate(differentiate(exp(t)),1) == exp(t)
    @test integrate(cos(t)) == sin(t)

    @test promote(ta(0.0), t) == (ta(0.0),ta(0.0))

    @test inverse(exp(t)-1) ≈ log(1+t)
    cfs = [(-n)^(n-1)/factorial(n) for n = 1:15]
    @test norm(inverse(t*exp(t))[1:end]./cfs .- 1) < 4tol1
    @test inverse(tan(t))(tan(t)) ≈ t
    @test atan(inverse(atan(t))) ≈ t
    @test inverse_map(sin(t))(sin(t)) ≈ t
    @test sinh(inverse_map(sinh(t))) ≈ t
    @test inverse_map(tanh(t)) ≈ inverse(tanh(t))

    @test_throws ArgumentError Taylor1([1,2,3], -2)
    @test_throws DomainError abs(ta(big(0)))
    @test_throws ArgumentError 1/t
    @test_throws ArgumentError zt/zt
    @test_throws DomainError t^1.5
    @test_throws ArgumentError t^(-2)
    @test_throws DomainError sqrt(t)
    @test_throws DomainError log(t)
    @test_throws ArgumentError cos(t)/sin(t)
    @test_throws AssertionError differentiate(30, exp(ta(1.0pi)))
    @test_throws DomainError inverse(exp(t))
    @test_throws DomainError abs(t)

    use_show_default(true)
    aa = sqrt(2)+Taylor1(2)
    @test string(aa) == "Taylor1{Float64}([1.4142135623730951, 1.0, 0.0], 2)"
    @test string([aa, aa]) ==
        "Taylor1{Float64}[Taylor1{Float64}([1.4142135623730951, 1.0, 0.0], 2), " *
        "Taylor1{Float64}([1.4142135623730951, 1.0, 0.0], 2)]"
    use_show_default(false)
    @test string(aa) == " 1.4142135623730951 + 1.0 t + 𝒪(t³)"
    set_taylor1_varname(" x ")
    @test string(aa) == " 1.4142135623730951 + 1.0 x + 𝒪(x³)"
    set_taylor1_varname(1, " t ")
    displayBigO(false)
    @test string(ta(-3)) == " - 3 + 1 t "
    @test string(ta(0)^3-3) == " - 3 + 1 t³ "
    @test TS.pretty_print(ta(3im)) == " ( 0 + 3im ) + ( 1 + 0im ) t "
    @test string(Taylor1([1,2,3,4,5], 2)) == string(Taylor1([1,2,3]))
    displayBigO(true)
    @test string(ta(-3)) == " - 3 + 1 t + 𝒪(t¹⁶)"
    @test string(ta(0)^3-3) == " - 3 + 1 t³ + 𝒪(t¹⁶)"
    @test TS.pretty_print(ta(3im)) == " ( 0 + 3im ) + ( 1 + 0im ) t + 𝒪(t¹⁶)"
    @test string(Taylor1([1,2,3,4,5], 2)) == string(Taylor1([1,2,3]))

    a = collect(1:12)
    t_a = Taylor1(a,15)
    t_C = complex(3.0,4.0) * t_a
    rnd = rand(10)
    @test typeof( norm(Taylor1(rnd)) ) == Float64
    @test norm(Taylor1(rnd)) > 0
    @test norm(t_a) == norm(a)
    @test norm(Taylor1(a,15), 3) == sum((a.^3))^(1/3)
    @test norm(t_a, Inf) == 12
    @test norm(t_C) == norm(complex(3.0,4.0)*a)

    @test TS.rtoldefault(Taylor1{Int}) == 0
    @test TS.rtoldefault(Taylor1{Float64}) == sqrt(eps(Float64))
    @test TS.rtoldefault(Taylor1{BigFloat}) == sqrt(eps(BigFloat))
    @test TS.real(Taylor1{Float64}) == Taylor1{Float64}
    @test TS.real(Taylor1{Complex{Float64}}) == Taylor1{Float64}
    @test isfinite(t_C)
    @test isfinite(t_a)
    @test !isfinite( Taylor1([0, Inf]) )
    @test !isfinite( Taylor1([NaN, 0]) )
    b = convert(Vector{Float64}, a)
    b[3] += eps(10.0)
    b[5] -= eps(10.0)
    t_b = Taylor1(b,15)
    t_C2 = t_C+eps(100.0)
    t_C3 = t_C+eps(100.0)*im
    @test isapprox(t_C, t_C)
    @test t_a ≈ t_a
    @test t_a ≈ t_b
    @test t_C ≈ t_C2
    @test t_C ≈ t_C3
    @test t_C3 ≈ t_C2
    t = Taylor1(25)
    p = sin(t)
    q = sin(t+eps())
    @test t ≈ t
    @test t ≈ t+sqrt(eps())
    @test isapprox(p, q, atol=eps())

    tf = Taylor1(35)
    @test Taylor1([180.0, rad2deg(1.0)], 35) == rad2deg(pi+tf)
    @test sin(pi/2+deg2rad(1.0)tf) == sin(deg2rad(90+tf))
    a = Taylor1(rand(10))
    b = Taylor1(rand(10))
    c = deepcopy(a)
    TS.deg2rad!(b, a, 0)
    @test a == c
    @test a[0]*(pi/180) == b[0]
    # TS.deg2rad!.(b, a, [0,1,2])
    # @test a == c
    # for i in 0:2
    #     @test a[i]*(pi/180) == b[i]
    # end
    a = Taylor1(rand(10))
    b = Taylor1(rand(10))
    c = deepcopy(a)
    TS.rad2deg!(b, a, 0)
    @test a == c
    @test a[0]*(180/pi) == b[0]
    # TS.rad2deg!.(b, a, [0,1,2])
    # @test a == c
    # for i in 0:2
    #     @test a[i]*(180/pi) == b[i]
    # end

    x = Taylor1([5.0,-1.5,3.0,-2.0,-20.0])
    @test x*x == x^2
    @test x*x*x == x*(x^2) == TaylorSeries.power_by_squaring(x, 3)
    @test x*x*x*x == (x^2)*(x^2) == TaylorSeries.power_by_squaring(x, 4)
    @test (x - 1.0)^2 == 1 - 2x + x^2
    @test (x - 1.0)^3 == -1 + 3x - 3x^2 + x^3
    @test (x - 1.0)^4 == 1 - 4x + 6x^2 - 4x^3 + x^4

    # Test additional Taylor1 constructors
    @test Taylor1{Float64}(true) == Taylor1([1.0])
    @test Taylor1{Float64}(false) == Taylor1([0.0])
    @test Taylor1{Int}(true) == Taylor1([1])
    @test Taylor1{Int}(false) == Taylor1([0])

    # Test fallback pretty_print
    st = Taylor1([SymbNumber(:x₀), SymbNumber(:x₁)])
    @test string(st) == " SymbNumber(:x₀) + SymbNumber(:x₁) t + 𝒪(t²)"

    @testset "Test Base.float overloads for Taylor1" begin
        @test float(Taylor1(-3, 2)) == Taylor1(-3.0, 2)
        @test float(Taylor1(-1//2, 2)) == Taylor1(-0.5, 2)
        @test float(Taylor1(3 - 0im, 2)) == Taylor1(3.0 - 0.0im, 2)
        x = Taylor1(rand(5))
        @test float(x) == x
        @test float(Taylor1{Int32}) == Taylor1{Float64}
        @test float(Taylor1{Int}) == Taylor1{Float64}
        @test float(Taylor1{Complex{Int}}) == Taylor1{ComplexF64}
    end
end

@testset "Test inv for Matrix{Taylor1{Float64}}" begin
    t = Taylor1(5)
    a = Diagonal(rand(0:10,3)) + rand(3, 3)
    ainv = inv(a)
    b = Taylor1.(a, 5)
    binv = inv(b)
    c = Symmetric(b)
    cinv = inv(c)
    tol = 1.0e-11

    for its = 1:10
        a .= Diagonal(rand(2:12,3)) + rand(3, 3)
        ainv .= inv(a)
        b .= Taylor1.(a, 5)
        binv .= inv(b)
        c .= Symmetric(Taylor1.(a, 5))
        cinv .= inv(c)
        @test norm(binv - ainv, Inf) ≤ tol
        @test norm(b*binv - I, Inf) ≤ tol
        @test norm(binv*b - I, Inf) ≤ tol
        @test norm(triu(b)*inv(UpperTriangular(b)) - I, Inf) ≤ tol
        @test norm(inv(LowerTriangular(b))*tril(b) - I, Inf) ≤ tol
        ainv .= inv(Symmetric(a))
        @test norm(cinv - ainv, Inf) ≤ tol
        @test norm(c*cinv - I, Inf) ≤ tol
        @test norm(cinv*c - I, Inf) ≤ tol
        @test norm(triu(c)*inv(UpperTriangular(c)) - I, Inf) ≤ tol
        @test norm(inv(LowerTriangular(c))*tril(c) - I, Inf) ≤ tol

        b .= b .+ t
        binv .= inv(b)
        @test norm(b*binv - I, Inf) ≤ tol
        @test norm(binv*b - I, Inf) ≤ tol
        @test norm(triu(b)*inv(triu(b)) - I, Inf) ≤ tol
        @test norm(inv(tril(b))*tril(b) - I, Inf) ≤ tol
        c .= Symmetric(b)
        cinv .= inv(c)
        @test norm(c*cinv - I, Inf) ≤ tol
        @test norm(cinv*c - I, Inf) ≤ tol
    end
end

@testset "Matrix multiplication for Taylor1" begin
    order = 30
    n1 = 100
    k1 = 90

    order = max(n1,k1)
    B1 = randn(n1,order)
    Y1 = randn(k1,order)

    A1  = randn(k1,n1)

    for A in (A1,sparse(A1))
        # B and Y contain elements of different orders
        B  = Taylor1{Float64}[Taylor1(collect(B1[i,1:i]),i) for i=1:n1]
        Y  = Taylor1{Float64}[Taylor1(collect(Y1[k,1:k]),k) for k=1:k1]
        Bcopy = deepcopy(B)
        mul!(Y,A,B)

        # do we get the same result when using the `A*B` form?
        @test A*B≈Y
        # Y should be extended after the multilpication
        @test reduce(&, [y1.order for y1 in Y] .== Y[1].order)
        # B should be unchanged
        @test B==Bcopy

        # is the result compatible with the matrix multiplication?  We
        # only check the zeroth order of the Taylor series.
        y1=sum(Y)[0]
        Y=A*B1[:,1]
        y2=sum(Y)

        # There is a small numerical error when comparing the generic
        # multiplication and the specialized version
        @test abs(y1-y2) < n1*(eps(y1)+eps(y2))

        @test_throws DimensionMismatch mul!(Y,A[:,1:end-1],B)
        @test_throws DimensionMismatch mul!(Y,A[1:end-1,:],B)
        @test_throws DimensionMismatch mul!(Y,A,B[1:end-1])
        @test_throws DimensionMismatch mul!(Y[1:end-1],A,B)
    end
end
