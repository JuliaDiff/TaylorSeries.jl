# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries
using Compat

if VERSION < v"0.7.0-DEV.2004"
    using Base.Test
    eeuler = Base.e
else
    using Test
    using LinearAlgebra, SparseArrays
    eeuler = Base.MathConstants.e
end

@testset "Tests for Taylor1 expansions" begin
    ta(a) = Taylor1([a,one(a)],15)
    t = Taylor1(Int,15)
    tim = im*t
    zt = zero(t)
    ot = 1.0*one(t)
    tol1 = eps(1.0)

    @test Taylor1 <: AbstractSeries
    @test Taylor1{Float64} <: AbstractSeries{Float64}

    @test Taylor1([1,2,3,4,5], 2) == Taylor1([1,2,3])
    @test get_order(Taylor1([1,2,3,4,5], 2)) == 2

    # ta = Taylor1([1,2,3,4,5])
    # vT = [h for h in ta]
    # voT = [h.order for h in Taylor1([1,2,3,4,5])]
    # @test vT == [Taylor1(1,0), Taylor1(2,1), Taylor1(3,2), Taylor1(4,3), Taylor1(5,4)]
    # @test voT == [0, 1, 2, 3, 4]

    v = [1,2]
    @test @compat typeof(TaylorSeries.resize_coeffs1!(v,3)) == Nothing
    @test v == [1,2,0,0]
    TaylorSeries.resize_coeffs1!(v,0)
    @test v == [1]
    TaylorSeries.resize_coeffs1!(v,3)
    setindex!(Taylor1(v),3,2)
    @test v == [1,0,3,0]
    @test Taylor1(v)[:] == [1,0,3,0]
    @test Taylor1(v)[:] == Taylor1(v).coeffs[:]
    setindex!(Taylor1(v),0,0:2)
    @test v == zero(v)
    setindex!(Taylor1(v),1,:)
    @test v == ones(Int, 4)
    Taylor1(v)[:] = 0
    @test v == zero(v)
    rv = [rand(0:3) for i in 1:4]
    @test Taylor1(rv)[:] == rv
    y = sin(Taylor1(16))
    @test y[:] == y.coeffs
    y[:] .= cos(Taylor1(16))[:]
    @test y == cos(Taylor1(16))
    @test y[:] == cos(Taylor1(16))[:]
    y = sin(Taylor1(16))
    rv = rand(5)
    y[0:4] = rv
    @test y[0:4] == rv
    @test y[5:end] == y.coeffs[6:end]
    rv = rand( length(y.coeffs) )
    y[:] = rv
    @test y[:] == rv
    y[:] = cos(Taylor1(16)).coeffs
    @test y == cos(Taylor1(16))
    @test y[:] == cos(Taylor1(16))[:]
    y[:] = 0.0
    @test y[:] == zero(y[:])
    y = sin(Taylor1(16))
    rv = rand.(length(0:4))
    y[0:4] .= rv
    @test y[0:4] == rv
    @test y[6:end] == sin(Taylor1(16))[6:end]
    rv = rand.(length(y))
    y[:] = rv
    @test y[:] == rv

    @test Taylor1([0,1,0,0]) == Taylor1(3)
    @test getcoeff(Taylor1(Complex{Float64},3),1) == complex(1.0,0.0)
    @test Taylor1(Complex{Float64},3)[1] == complex(1.0,0.0)
    @test getindex(Taylor1(3),1) == 1.0
    @inferred convert(Taylor1{Complex{Float64}},ot) == Taylor1{Complex{Float64}}
    @test eltype(convert(Taylor1{Complex{Float64}},ot)) == Complex{Float64}
    @test eltype(convert(Taylor1{Complex{Float64}},1)) == Complex{Float64}
    @test eltype(convert(Taylor1, 1im)) == Complex{Int}
    @test convert(Taylor1, 1im) == Taylor1(1im, 0)
    @test convert(Taylor1{Int},[0,2]) == 2*t
    @test convert(Taylor1{Complex{Int}},[0,2]) == (2+0im)*t
    @test convert(Taylor1{BigFloat},[0.0, 1.0]) == ta(big(0.0))
    @test promote(t,Taylor1(1.0,0)) == (ta(0.0),ot)
    @test promote(0,Taylor1(1.0,0)) == (zt,ot)
    @test eltype(promote(ta(0),zeros(Int,2))[2]) == Int
    @test eltype(promote(ta(0.0),zeros(Int,2))[2]) == Float64
    @test eltype(promote(0,Taylor1(ot))[1]) == Float64
    @test eltype(promote(1.0+im, zt)[1]) == Complex{Float64}

    @test length(Taylor1(10)) == 11
    @test length(TaylorSeries.fixorder(zt,Taylor1([1]))[2]) == 16
    @test eltype(TaylorSeries.fixorder(zt,Taylor1([1]))[1]) == Int
    @test TaylorSeries.findfirst(t) == 1
    @test TaylorSeries.findfirst(zt) == -1
    @test iszero(zero(t))
    @test !iszero(one(t))
    @test isinf(Taylor1([typemax(1.0)]))
    @test isnan(Taylor1([typemax(1.0), NaN]))

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
    @test one(t)/one(t) == 1.0
    @test tsquare/t == t
    @test t/(t*3) == (1/3)*ot
    @test t/3im == -tim/3
    @test 1/(1-t) == Taylor1(ones(t.order+1))
    @test Taylor1([0,1,1])/t == t+1
    @test (t+im)^2 == tsquare+2im*t-1
    @test (t+im)^3 == Taylor1([-1im,-3,3im,1],15)
    @test (t+im)^4 == Taylor1([1,-4im,-6,4im,1],15)
    @test imag(tsquare+2im*t-1) == 2t
    @test (Rational(1,2)*tsquare)[2] == 1//2
    @test t^2/tsquare == ot
    @test ((1+t)^(1/3))[2]+1/9 ≤ tol1
    @test 1-tsquare == (1+t)-t*(1+t)
    @test (1-tsquare)^2 == (1+t)^2.0 * (1-t)^2.0
    @test (sqrt(1+t))[2] == -1/8
    @test ((1-tsquare)^(1//2))^2 == 1-tsquare
    @test ((1-t)^(1//4))[14] == -4188908511//549755813888
    @test abs(((1+t)^3.2)[13] + 5.4021062656e-5) < tol1
    @test Taylor1(BigFloat,5)/6 == 1im*Taylor1(5)/complex(0,BigInt(6))
    @test Taylor1(BigFloat,5)/(6*Taylor1(3)) == 1/BigInt(6)
    @test Taylor1(BigFloat,5)/(6im*Taylor1(3)) == -1im/BigInt(6)

    # These tests involve some sort of factorization
    @test t/(t+t^2) == 1/(1+t)
    @test sqrt(t^2+t^3) == t*sqrt(1+t)
    @test (t^3+t^4)^(1/3) ≈ t*(1+t)^(1/3)
    @test norm((t^3+t^4)^(1/3) - t*(1+t)^(1/3), Inf) < eps()
    @test ((t^3+t^4)^(1/3))[15] ≈ -8617640/1162261467

    trational = ta(0//1)
    @inferred ta(0//1) == Taylor1{Rational{Int}}
    @test eltype(trational) == Rational{Int}
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
    @test eltype(taylor_expand(x->x^2+1,1)) == eltype(1)
    tsq = t^2
    update!(tsq,2.0)
    @test tsq == (t+2.0)^2
    update!(tsq,-2.0)
    @test tsq == t^2

    @test log(exp(tsquare)) == tsquare
    @test exp(log(1-tsquare)) == 1-tsquare
    @test log((1-t)^2) == 2*log(1-t)
    @test real(exp(tim)) == cos(t)
    @test imag(exp(tim)) == sin(t)
    @test exp(conj(tim)) == cos(t)-im*sin(t) == exp(tim')
    @test (exp(t))^(2im) == cos(2t)+im*sin(2t)
    @test (exp(t))^Taylor1([-5.2im]) == cos(5.2t)-im*sin(5.2t)
    @test getcoeff(convert(Taylor1{Rational{Int}},cos(t)),8) == 1//factorial(8)
    @test abs((tan(t))[7]- 17/315) < tol1
    @test abs((tan(t))[13]- 21844/6081075) < tol1
    @test evaluate(exp(Taylor1([0,1],17)),1.0) == 1.0*eeuler
    @test evaluate(exp(Taylor1([0,1],1))) == 1.0
    @test evaluate(exp(t),t^2) == exp(t^2)
    #Test function-like behavior for Taylor1s
    t17 = Taylor1([0,1],17)
    myexpfun = exp(t17)
    @test myexpfun(1.0) == 1.0*eeuler
    @test myexpfun() == 1.0
    @test myexpfun(t17^2) == exp(t17^2)
    @test exp(t17^2)(t17) == exp(t17^2)
    p = cos(t17)
    q = sin(t17)
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
    vr = rand(2)
    @test p.(vr) == evaluate.(p, vr)
    Mr = rand(3,3,3)
    @test p.(Mr) == evaluate.(p, Mr)
    mytaylor1 = Taylor1(rand(20))
    vr = rand(5)
    @test p(vr) == p.(vr)
    @test p(vr) == evaluate.(p,vr)
    @test p(Mr) == p.(Mr)
    @test p(Mr) == evaluate.(p,Mr)
    taylor_a = Taylor1(Int64,10)
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
    @test asin(t) + acos(t) == pi/2
    @test derivative(acos(t)) == - 1/sqrt(1-t^2)

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
    TaylorSeries.add!(tt, ut, ut, 1)
    @test tt[1] == 2.0
    TaylorSeries.add!(tt, -3.0, 0)
    @test tt[0] == -3.0
    TaylorSeries.add!(tt, -3.0, 1)
    @test tt[1] == 0.0
    TaylorSeries.subst!(tt, ut, ut, 1)
    @test tt[1] == 0.0
    TaylorSeries.subst!(tt, -3.0, 0)
    @test tt[0] == 3.0
    TaylorSeries.subst!(tt, -2.5, 1)
    @test tt[1] == 0.0
    iind, cind = TaylorSeries.divfactorization(ut, ut)
    @test iind == 1
    @test cind == 1.0
    TaylorSeries.div!(tt, ut, ut, 0, iind)
    @test tt[0] == cind
    TaylorSeries.div!(tt, 1+ut, 1+ut, 0)
    @test tt[0] == 1.0
    TaylorSeries.pow!(tt, 1+t, 1.5, 0, 0)
    @test tt[0] == 1.0
    TaylorSeries.pow!(tt, 1+t, 1.5, 0)
    @test tt[0] == 1.0
    TaylorSeries.sqrt!(tt, 1+t, 0, 0)
    @test tt[0] == 1.0
    TaylorSeries.sqrt!(tt, 1+t, 0)
    @test tt[0] == 1.0
    TaylorSeries.exp!(tt, t, 0)
    @test tt[0] == exp(t[0])
    TaylorSeries.log!(tt, 1.0+t, 0)
    @test tt[0] == 0.0
    ct = zero(ut)
    TaylorSeries.sincos!(tt, ct, t, 0)
    @test tt[0] == sin(t[0])
    @test ct[0] == cos(t[0])
    TaylorSeries.tan!(tt, t, ct, 0)
    @test tt[0] == tan(t[0])
    @test ct[0] == tan(t[0])^2
    TaylorSeries.asin!(tt, t, ct, 0)
    @test tt[0] == asin(t[0])
    @test ct[0] == sqrt(1.0-t[0]^2)
    TaylorSeries.acos!(tt, t, ct, 0)
    @test tt[0] == acos(t[0])
    @test ct[0] == sqrt(1.0-t[0]^2)
    TaylorSeries.atan!(tt, ut, ct, 0)
    @test tt[0] == atan(t[0])
    @test ct[0] == 1.0+t[0]^2
    TaylorSeries.sinhcosh!(tt, ct, ut, 0)
    @test tt[0] == sinh(t[0])
    @test ct[0] == cosh(t[0])
    TaylorSeries.tanh!(tt, ut, ct, 0)
    @test tt[0] == tanh(t[0])
    @test ct[0] == tanh(t[0])^2

    v = [sin(t), exp(-t)]
    @compat vv = Vector{Float64}(uninitialized, 2)
    @test evaluate!(v, zero(Int), vv) == nothing
    @test vv == [0.0,1.0]
    @test evaluate(v) == vv
    @test evaluate(v, complex(0.0,0.2)) ==
        [complex(0.0,sinh(0.2)),complex(cos(0.2),sin(-0.2))]

    @test derivative(exp(ta(1.0)), 0) == exp(ta(1.0))
    expected_result_approx = Taylor1(convert(Vector{Float64},exp(ta(1.0))[0:10]))
    @test derivative(exp(ta(1.0)), 5) ≈ expected_result_approx atol=eps() rtol=0.0
    expected_result_approx = Taylor1(convert(Vector{Float64},exp(ta(1.0pi))[0:12]),15)
    @test derivative(exp(ta(1.0pi)), 3) ≈ expected_result_approx atol=eps(16.0) rtol=0.0
    expected_result_approx = Taylor1(convert(Vector{Float64},exp(ta(1.0pi))[0:5]),15)
    @test derivative(exp(ta(1.0pi)), 10) ≈ expected_result_approx atol=eps(64.0) rtol=0.0



    @test derivative(exp(ta(1.0)), 5)() == exp(1.0)
    @test derivative(exp(ta(1.0pi)), 3)() == exp(1.0pi)
    @test isapprox(derivative(exp(ta(1.0pi)), 10)() , exp(1.0pi) )

    @test derivative(5, exp(ta(1.0))) == exp(1.0)
    @test derivative(3, exp(ta(1.0pi))) == exp(1.0pi)
    @test isapprox(derivative(10, exp(ta(1.0pi))) , exp(1.0pi) )

    @test integrate(derivative(exp(t)),1) == exp(t)
    @test integrate(cos(t)) == sin(t)

    @test promote(ta(0.0), t) == (ta(0.0),ta(0.0))

    @test norm((inverse(exp(t)-1) - log(1+t)).coeffs) < 2tol1
    cfs = [(-n)^(n-1)/gamma(n+1) for n = 1:15]
    @test norm(inverse(t*exp(t))[1:end]./cfs .- 1) < 4tol1

    @test_throws ArgumentError Taylor1([1,2,3], -2)
    @test_throws ArgumentError abs(ta(big(0)))
    @test_throws ArgumentError 1/t
    @test_throws ArgumentError zt/zt
    @test_throws ArgumentError t^1.5
    if VERSION < v"0.7.0-DEV"
        @test_throws DomainError t^(-2)
    else
        @test_throws ArgumentError t^(-2)
    end
    @test_throws ArgumentError sqrt(t)
    @test_throws ArgumentError log(t)
    @test_throws ArgumentError cos(t)/sin(t)
    @test_throws AssertionError derivative(30, exp(ta(1.0pi)))
    @test_throws ArgumentError inverse(exp(t))

    displayBigO(false)
    @test string(ta(-3)) == " - 3 + 1 t "
    @test string(ta(0)^3-3) == " - 3 + 1 t³ "
    @test TaylorSeries.pretty_print(ta(3im)) == " ( 3 im )  + ( 1 ) t "
    @test string(Taylor1([1,2,3,4,5], 2)) == string(Taylor1([1,2,3]))
    displayBigO(true)
    @test string(ta(-3)) == " - 3 + 1 t + 𝒪(t¹⁶)"
    @test string(ta(0)^3-3) == " - 3 + 1 t³ + 𝒪(t¹⁶)"
    @test TaylorSeries.pretty_print(ta(3im)) == " ( 3 im )  + ( 1 ) t + 𝒪(t¹⁶)"
    @test string(Taylor1([1,2,3,4,5], 2)) == string(Taylor1([1,2,3]))

    a = collect(1:12)
    t_a = Taylor1(a,15)
    t_C = complex(3.0,4.0) * t_a
    rnd = rand(10)
    @test typeof( norm(Taylor1(rnd)) ) == Float64
    @test norm(Taylor1(rnd)) > 0
    @test norm(t_a) == norm(a)
    @test norm(Taylor1(a,15),3) == sum((a.^3))^(1/3)
    @test norm(t_a,Inf) == 12
    @test norm(t_C) == norm(complex(3.0,4.0)*a)

    @test TaylorSeries.rtoldefault(Taylor1{Int64}) == 0
    @test TaylorSeries.rtoldefault(Taylor1{Float64}) == sqrt(eps(Float64))
    @test TaylorSeries.rtoldefault(Taylor1{BigFloat}) == sqrt(eps(BigFloat))
    @test TaylorSeries.real(Taylor1{Float64}) == Taylor1{Float64}
    @test TaylorSeries.real(Taylor1{Complex{Float64}}) == Taylor1{Float64}
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

    t = Taylor1(35)
    @test Taylor1([180.0, rad2deg(1.0)], 35) == rad2deg(pi+t)
    @test sin(pi/2+deg2rad(1.0)t) == sin(deg2rad(90+t))
    a = Taylor1(rand(10))
    b = Taylor1(rand(10))
    c = deepcopy(a)
    TaylorSeries.deg2rad!(b, a, 0)
    @test a == c
    @test a[0]*(pi/180) == b[0]
    TaylorSeries.deg2rad!.(b, a, [0,1,2])
    @test a == c
    for i in 0:2
        @test a[i]*(pi/180) == b[i]
    end
    a = Taylor1(rand(10))
    b = Taylor1(rand(10))
    c = deepcopy(a)
    TaylorSeries.rad2deg!(b, a, 0)
    @test a == c
    @test a[0]*(180/pi) == b[0]
    TaylorSeries.rad2deg!.(b, a, [0,1,2])
    @test a == c
    for i in 0:2
        @test a[i]*(180/pi) == b[i]
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
        @compat mul!(Y,A,B)

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

        @compat @test_throws DimensionMismatch mul!(Y,A[:,1:end-1],B)
        @compat @test_throws DimensionMismatch mul!(Y,A[1:end-1,:],B)
        @compat @test_throws DimensionMismatch mul!(Y,A,B[1:end-1])
        @compat @test_throws DimensionMismatch mul!(Y[1:end-1],A,B)
    end
end
