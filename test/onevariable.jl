# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries
using Base.Test

@testset "Tests for Taylor1 expansions" begin
    ta(a) = Taylor1([a,one(a)],15)
    t = Taylor1(Int,15)
    tim = im*t
    zt = zero(t)
    ot = 1.0*one(t)
    tol1 = eps(1.0)

    @test Taylor1 <: AbstractSeries
    @test Taylor1{Float64} <: AbstractSeries{Float64}

    @test Taylor1([1,2,3,4,5], 2) == Taylor1([1,2,3,4,5])

    v = [1,2]
    @test typeof(TaylorSeries.resize_coeffs1!(v,3)) == Void
    @test v == [1,2,0,0]
    TaylorSeries.resize_coeffs1!(v,0)
    @test v == [1,2,0,0]
    setindex!(Taylor1(v),3,3)
    @test v == [1,2,3,0]
    @test Taylor1(v)[:] == [1,2,3,0]
    @test Taylor1(v)[:] == Taylor1(v).coeffs[:]
    setindex!(Taylor1(v),0,1:3)
    @test v == zero(v)
    setindex!(Taylor1(v),1,:)
    @test v == [1,1,1,1]
    Taylor1(v)[:] = 0
    @test v == zero(v)
    rv = [rand(0:3) for i in 1:4]
    @test Taylor1(rv)[:] == Taylor1(rv)[1:end]
    y = sin(Taylor1(16))
    @test y[:] == sin(Taylor1(16))[:]
    y[:] .= cos(Taylor1(16))[:]
    @test y == cos(Taylor1(16))
    @test y[:] == cos(Taylor1(16))[:]
    y = sin(Taylor1(16))
    rv = rand(5)
    y[1:5] = rv
    @test y[1:5] == rv
    @test y[6:end] == sin(Taylor1(16))[6:end]
    rv = rand( length(y.coeffs) )
    y[:] = rv
    @test y[:] == rv
    y[:] = cos(Taylor1(16)).coeffs
    @test y == cos(Taylor1(16))
    @test y[:] == cos(Taylor1(16))[:]
    y[:] = 0.0
    @test y[:] == zero(y[:])
    y = sin(Taylor1(16))
    rv = y[1:5] .= rand.()
    @test y[1:5] == rv
    @test y[6:end] == sin(Taylor1(16))[6:end]
    rv = y[:] .= rand.()
    @test y[:] == rv

    @test Taylor1([0,1,0,0]) == Taylor1(3)
    @test get_coeff(Taylor1(Complex128,3),1) == complex(1.0,0.0)
    @test Taylor1(Complex128,3)[2] == complex(1.0,0.0)
    @test getindex(Taylor1(3),2) == 1.0
    @inferred convert(Taylor1{Complex128},ot) == Taylor1{Complex128}
    @test eltype(convert(Taylor1{Complex128},ot)) == Complex128
    @test eltype(convert(Taylor1{Complex128},1)) == Complex128
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
    @test get_coeff(tim,1) == complex(0,1)
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
    @test (Rational(1,2)*tsquare)[3] == 1//2
    @test t^2/tsquare == ot
    @test ((1+t)^(1/3))[3]+1/9 ≤ tol1
    @test 1-tsquare == (1+t)-t*(1+t)
    @test (1-tsquare)^2 == (1+t)^2.0 * (1-t)^2.0
    @test (sqrt(1+t))[3] == -1/8
    @test ((1-tsquare)^(1//2))^2 == 1-tsquare
    @test ((1-t)^(1//4))[15] == -4188908511//549755813888
    @test abs(((1+t)^3.2)[14] + 5.4021062656e-5) < tol1
    @test Taylor1(BigFloat,5)/6 == 1im*Taylor1(5)/complex(0,BigInt(6))
    @test Taylor1(BigFloat,5)/(6*Taylor1(3)) == 1/BigInt(6)
    @test Taylor1(BigFloat,5)/(6im*Taylor1(3)) == -1im/BigInt(6)

    trational = ta(0//1)
    @inferred ta(0//1) == Taylor1{Rational{Int}}
    @test eltype(trational) == Rational{Int}
    @test trational + 1//3 == Taylor1([1//3,1],15)
    @test complex(3,1)*trational^2 == Taylor1([0//1,0//1,complex(3,1)//1],15)
    @test trational^2/3 == Taylor1([0//1,0//1,1//3],15)
    @test trational^3/complex(7,1) == Taylor1([0,0,0,complex(7//50,-1//50)],15)
    @test sqrt(zero(t)) == zero(t)

    @test isapprox( rem(4.1 + t,4)[1], 0.1 )
    @test isapprox( mod(4.1 + t,4)[1], 0.1 )
    @test isapprox( rem(1+Taylor1(Int,4),4.0)[1], 1.0 )
    @test isapprox( mod(1+Taylor1(Int,4),4.0)[1], 1.0 )
    @test isapprox( mod2pi(2pi+0.1+t)[1], 0.1 )

    @test abs(ta(1)) == ta(1)
    @test abs(ta(-1.0)) == -ta(-1.0)

    @test log(exp(tsquare)) == tsquare
    @test exp(log(1-tsquare)) == 1-tsquare
    @test log((1-t)^2) == 2*log(1-t)
    @test real(exp(tim)) == cos(t)
    @test imag(exp(tim)) == sin(t)
    @test exp(conj(tim)) == cos(t)-im*sin(t) == exp(tim')
    @test (exp(t))^(2im) == cos(2t)+im*sin(2t)
    @test (exp(t))^Taylor1([-5.2im]) == cos(5.2t)-im*sin(5.2t)
    @test get_coeff(convert(Taylor1{Rational{Int}},cos(t)),8) == 1//factorial(8)
    @test abs((tan(t))[8]- 17/315) < tol1
    @test abs((tan(t))[14]- 21844/6081075) < tol1
    @test evaluate(exp(Taylor1([0,1],17)),1.0) == 1.0*e
    @test evaluate(exp(Taylor1([0,1],1))) == 1.0
    @test evaluate(exp(t),t^2) == exp(t^2)

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
    @test tt[2] == 2.0
    TaylorSeries.subst!(tt, ut, ut, 1)
    @test tt[2] == 0.0
    iind, cind = TaylorSeries.divfactorization(ut, ut)
    @test iind == 1
    @test cind == 1.0
    TaylorSeries.div!(tt, ut, ut, 0, iind)
    @test tt[1] == cind
    TaylorSeries.div!(tt, 1+ut, 1+ut, 0)
    @test tt[1] == 1.0
    TaylorSeries.pow!(tt, 1+t, 1.5, 0, 0)
    @test tt[1] == 1.0
    TaylorSeries.pow!(tt, 1+t, 1.5, 0)
    @test tt[1] == 1.0
    TaylorSeries.sqrt!(tt, 1+t, 0, 0)
    @test tt[1] == 1.0
    TaylorSeries.sqrt!(tt, 1+t, 0)
    @test tt[1] == 1.0
    TaylorSeries.exp!(tt, t, 0)
    @test tt[1] == exp(t[1])
    TaylorSeries.log!(tt, 1.0+t, 0)
    @test tt[1] == 0.0
    ct = zero(ut)
    TaylorSeries.sincos!(tt, ct, t, 0)
    @test tt[1] == sin(t[1])
    @test ct[1] == cos(t[1])
    TaylorSeries.tan!(tt, t, ct, 0)
    @test tt[1] == tan(t[1])
    @test ct[1] == tan(t[1])^2
    TaylorSeries.asin!(tt, t, ct, 0)
    @test tt[1] == asin(t[1])
    @test ct[1] == sqrt(1.0-t[1]^2)
    TaylorSeries.acos!(tt, t, ct, 0)
    @test tt[1] == acos(t[1])
    @test ct[1] == sqrt(1.0-t[1]^2)
    TaylorSeries.atan!(tt, ut, ct, 0)
    @test tt[1] == atan(t[1])
    @test ct[1] == 1.0+t[1]^2
    TaylorSeries.sinhcosh!(tt, ct, ut, 0)
    @test tt[1] == sinh(t[1])
    @test ct[1] == cosh(t[1])
    TaylorSeries.tanh!(tt, ut, ct, 0)
    @test tt[1] == tanh(t[1])
    @test ct[1] == tanh(t[1])^2

    v = [sin(t), exp(-t)]
    vv = Vector{Float64}(2)
    @test evaluate!(v, zero(Int), vv) == nothing
    @test vv == [0.0,1.0]
    @test evaluate(v) == vv
    @test evaluate(v, complex(0.0,0.2)) ==
        [complex(0.0,sinh(0.2)),complex(cos(0.2),sin(-0.2))]
    @test derivative(5, exp(ta(1.0))) == exp(1.0)
    @test derivative(3, exp(ta(1.0pi))) == exp(1.0pi)
    @test isapprox(derivative(10, exp(ta(1.0pi))) , exp(1.0pi) )
    @test integrate(derivative(exp(t)),1) == exp(t)
    @test integrate(cos(t)) == sin(t)

    @test promote(ta(0.0), t) == (ta(0.0),ta(0.0))

    @test norm((inverse(exp(t)-1) - log(1+t)).coeffs) < 2tol1
    cfs = [(-n)^(n-1)/gamma(n+1) for n = 1:15]
    @test norm(inverse(t*exp(t))[2:end]./cfs-1) < 4tol1

    @test_throws ArgumentError abs(ta(big(0)))
    @test_throws ArgumentError 1/t
    @test_throws ArgumentError zt/zt
    @test_throws ArgumentError t^1.5
    @test_throws DomainError t^(-2)
    @test_throws ArgumentError sqrt(t)
    @test_throws ArgumentError log(t)
    @test_throws ArgumentError cos(t)/sin(t)
    @test_throws AssertionError derivative(30, exp(ta(1.0pi)))
    @test_throws ArgumentError inverse(exp(t))

    @test string(ta(-3)) == " - 3 + 1 t + 𝒪(t¹⁶)"
    @test string(ta(0)^3-3) == " - 3 + 1 t³ + 𝒪(t¹⁶)"
    @test TaylorSeries.pretty_print(ta(3im)) == " ( 3 im )  + ( 1 ) t + 𝒪(t¹⁶)"
    @test string(Taylor1([1,2,3,4,5], 2)) == string(Taylor1([1,2,3,4,5]))
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
        A_mul_B!(Y,A,B)

        # do we get the same result when using the `A*B` form?
        @test A*B==Y
        # Y should be extended after the multilpication
        @test reduce(&, [y1.order for y1 in Y] .== Y[1].order)
        # B should be unchanged
        @test B==Bcopy

        # is the result compatible with the matrix multiplication?  We
        # only check the zeroth order of the Taylor series.
        y1=sum(Y)[1]
        Y=A*B1[:,1]
        y2=sum(Y)

        # There is a small numerical error when comparing the generic
        # multiplication and the specialized version
        @test abs(y1-y2) < n1*(eps(y1)+eps(y2))

        @test_throws DimensionMismatch A_mul_B!(Y,A[:,1:end-1],B)
        @test_throws DimensionMismatch A_mul_B!(Y,A[1:end-1,:],B)
        @test_throws DimensionMismatch A_mul_B!(Y,A,B[1:end-1])
        @test_throws DimensionMismatch A_mul_B!(Y[1:end-1],A,B)
    end
end
