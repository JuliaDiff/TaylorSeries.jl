# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries
using Base.Test

@testset "Tests for HomogeneousPolynomial and TaylorN" begin
    @test HomogeneousPolynomial <: AbstractSeries
    @test HomogeneousPolynomial{Int} <: AbstractSeries{Int}
    @test TaylorN{Float64} <: AbstractSeries{Float64}

    @test eltype(set_variables(Int, "x", numvars=2, order=6)) == TaylorN{Int}
    @test eltype(set_variables("x", numvars=2, order=6)) == TaylorN{Float64}
    @test eltype(set_variables(BigInt, "x y", order=6)) == TaylorN{BigInt}
    @test eltype(set_variables("x y", order=6)) == TaylorN{Float64}
    @test typeof(show_params_TaylorN()) == Void

    @test TaylorSeries.coeff_table[2][1] == [1,0]
    @test TaylorSeries.index_table[2][1] == 7
    @test TaylorSeries.in_base(get_order(),[2,1]) == 15
    @test TaylorSeries.pos_table[4][15] == 2

    @test get_order() == 6
    @test get_numvars() == 2

    x, y = set_variables("x y", order=6)
    @test x.order == 6
    @test TaylorSeries.set_variable_names(["x","y"]) == ["x", "y"]
    @test TaylorSeries.get_variable_names() == ["x", "y"]
    @test x == HomogeneousPolynomial(Float64, 1)
    @test x == HomogeneousPolynomial(1)
    @test y == HomogeneousPolynomial(Float64, 2)
    @test y == HomogeneousPolynomial(2)
    @test !isnan(x)

    set_variables("x", numvars=2, order=17)
    v = [1,2]
    @test typeof(TaylorSeries.resize_coeffsHP!(v,2)) == Void
    @test v == [1,2,0]
    @test_throws AssertionError TaylorSeries.resize_coeffsHP!(v,1)
    HomogeneousPolynomial(v)[3] = 3
    @test v == [1,2,3]
    HomogeneousPolynomial(v)[1:3] = 3
    @test v == [3,3,3]

    xH = HomogeneousPolynomial([1,0])
    yH = HomogeneousPolynomial([0,1],1)
    @test HomogeneousPolynomial(0,0)  == 0
    xT = TaylorN(xH, 17)
    yT = TaylorN(Int64, 2, order=17)
    zeroT = zero( TaylorN([xH],1) )
    @test zeroT.coeffs == zeros(HomogeneousPolynomial{Int}, 1)
    @test length(zeros(HomogeneousPolynomial{Int}, 1)) == 2
    @test one(HomogeneousPolynomial(1,1)) == HomogeneousPolynomial([1,1])
    uT = one(convert(TaylorN{Float64},yT))
    @test uT == one(HomogeneousPolynomial)
    @test zeroT[1] == HomogeneousPolynomial(0, 0)
    @test uT[1] == HomogeneousPolynomial(1, 0)
    @test ones(xH,1) == [1, xH+yH]
    @test typeof(ones(xH,2)) == Array{HomogeneousPolynomial{Int},1}
    @test length(ones(xH,2)) == 3
    @test ones(HomogeneousPolynomial{Complex{Int}},0) ==
        [HomogeneousPolynomial([complex(1,0)], 0)]
    @test !isnan(uT)

    @test get_order(zeroT) == 1
    @test xT[2][1] == 1
    @test yH[2] == 1
    @test get_coeff(xT,[1,0]) == 1
    @test get_coeff(yH,[1,0]) == 0
    @test typeof(convert(HomogeneousPolynomial,1im)) ==
        HomogeneousPolynomial{Complex{Int}}
    @test convert(HomogeneousPolynomial,1im) ==
        HomogeneousPolynomial([complex(0,1)], 0)
    @test convert(HomogeneousPolynomial{Int64},[1,1]) == xH+yH
    @test convert(HomogeneousPolynomial{Float64},[2,-1]) == 2.0xH-yH
    @test typeof(convert(TaylorN,1im)) == TaylorN{Complex{Int}}
    @test convert(TaylorN, 1im) ==
        TaylorN([HomogeneousPolynomial([complex(0,1)], 0)], 0)
    @test convert(TaylorN{Float64}, yH) == 1.0*yT
    @test convert(TaylorN{Float64}, [xH,yH]) == xT+1.0*yT
    @test convert(TaylorN{Int}, [xH,yH]) == xT+yT
    @test promote(xH, [1,1])[2] == xH+yH
    @test promote(xH, yT)[1] == xT
    @test promote(xT, [xH,yH])[2] == xT+yT
    @test typeof(promote(im*xT,[xH,yH])[2]) == TaylorN{Complex{Int64}}
    @test iszero(zeroT.coeffs)
    @test iszero(zero(xH))
    @test !iszero(uT)
    @test iszero(zeroT)

    @test HomogeneousPolynomial(xH,1) == HomogeneousPolynomial(xH)
    @test eltype(xH) == Int
    @test length(xH) == 2
    @test zero(xH) == 0*xH
    @test one(yH) == xH+yH
    @test xH * true == xH
    @test false * yH == zero(yH)
    @test get_order(yH) == 1
    @test get_order(xT) == 17
    @test xT * true == xT
    @test false * yT == zero(yT)

    @test xT == TaylorN([xH])
    @test one(xT) == TaylorN(1,5)
    @test TaylorN(zeroT,5) == 0
    @test TaylorN(uT) == convert(TaylorN{Complex},1)
    @test get_numvars() == 2
    @test length(uT) == get_order()+1
    @test eltype(convert(TaylorN{Complex128},1)) == Complex128

    @test 1+xT+yT == TaylorN(1,1) + TaylorN([xH,yH],1)
    @test xT-yT-1 == TaylorN([-1,xH-yH])
    @test xT*yT == TaylorN([HomogeneousPolynomial([0,1,0],2)])
    @test (1/(1-xT))[4] == HomogeneousPolynomial([1.0],3)
    @test xH^20 == HomogeneousPolynomial([0], get_order())
    @test (yT/(1-xT))[5] == xH^3 * yH
    @test mod(1+xT,1) == +xT
    @test (rem(1+xT,1))[1] == 0
    @test mod(1+xT,1.0) == +xT
    @test (rem(1+xT,1.0))[1] == 0
    @test abs(1-xT)  == 1-xT
    @test abs(-1-xT)  == 1+xT
    @test derivative(mod2pi(2pi+yT^3),2) == derivative(yT^3,2)
    @test derivative(yT) == zeroT
    @test -xT/3im == im*xT/3
    @test (xH/3im)' == im*xH/3
    @test xT/BigInt(3) == TaylorN(BigFloat,1)/3
    @test xT/complex(0,BigInt(3)) == -im*xT/BigInt(3)
    @test (xH/complex(0,BigInt(3)))' ==
        im*HomogeneousPolynomial([BigInt(1),0])/3
    @test evaluate(xH) == zero(eltype(xH))

    @test derivative(2xT*yT^2,1) == 2yT^2
    @test xT*xT^3 == xT^4
    txy = 1.0 + xT*yT - 0.5*xT^2*yT + (1/3)*xT^3*yT + 0.5*xT^2*yT^2
    @test getindex((1+TaylorN(1))^TaylorN(2),1:5) == txy.coeffs[1:5]
    @test ( (1+TaylorN(1))^TaylorN(2) )[:] == ( (1+TaylorN(1))^TaylorN(2) ).coeffs[:]
    @test txy.coeffs[:] == txy[:]
    @test txy.coeffs[:] == txy[1:end]
    txy[:] .= ( -1.0 + 3xT*yT - xT^2*yT + (4/3)*xT^3*yT + (1/3)*xT*yT^3 + 0.5*xT^2*yT^2 + 0.5*xT*yT^3 )[:]
    @test txy[:] == ( -1.0 + 3xT*yT - xT^2*yT + (4/3)*xT^3*yT + (1/3)*xT*yT^3 + 0.5*xT^2*yT^2 + 0.5*xT*yT^3 )[:]
    txy[2:end-1] .= ( 1.0 - xT*yT + 0.5*xT^2*yT - (2/3)*xT*yT^3 - 0.5*xT^2*yT^2  + 7*xT^3*yT )[2:end-1]
    @test txy[2:end-1] == ( 1.0 - xT*yT + 0.5*xT^2*yT - (2/3)*xT*yT^3 - 0.5*xT^2*yT^2  + 7*xT^3*yT )[2:end-1]

    a = -5.0 + sin(xT+yT^2)
    b = deepcopy(a)
    @test a[:] == a[1:end]
    @test a[:] == b[:]
    @test a[1:end] == b[1:end]
    @test a[end][:] == b[end][:]
    @test a[end][1:end] == b[end][1:end]
    rv = a[end][:] .= rand.()
    @test a[end][:] == rv
    @test a[end][:] != b[end][:]
    rv = a[end][1:end] .= rand.()
    @test a[end][1:end] == rv
    @test a[end][1:end] != b[end][1:end]

    hp = HomogeneousPolynomial(1)^8
    rv1 = rand( length(hp) )
    hp[:] = rv1
    @test rv1 == hp[:]
    rv2 = rand( length(hp)-2 )
    hp[1:end-2] = rv2
    @test hp[1:end-2] == rv2
    @test hp[end-1:end] == rv1[end-1:end]
    hp[3:4] = 0.0
    @test hp[1:2] == rv2[1:2]
    @test hp[3:4] == zeros(2)
    @test hp[5:end-2] == rv2[5:end]
    @test hp[end-1:end] == rv1[end-1:end]
    hp[:] = 0.0
    @test hp[:] == zero(rv1)

    pol = sin(xT+yT*xT)+yT^2-(1-xT)^3
    q = deepcopy(pol)
    q[:] = 0.0
    @test q[:] == zero(q[:])
    q[:] = pol.coeffs
    @test q == pol
    @test q[:] == pol[:]
    q[2:end-1] = 0.0
    @test q[2:end-1] == zero(q[2:end-1])
    @test q[1] == pol[1]
    @test q[end] == pol[end]
    q[:] = pol.coeffs
    @test q == pol
    @test q[:] == pol[:]
    zH0 = zero(HomogeneousPolynomial{Float64})
    q[:] = zH0
    @test q[:] == zero(q[:])
    q[:] = pol.coeffs
    @test q == pol
    @test q[:] == pol[:]
    q[2:end-1] = zH0
    @test q[2:end-1] == zero(q[2:end-1])
    @test q[1] == pol[1]
    @test q[end] == pol[end]


    @test_throws DomainError yT^(-2)
    @test_throws DomainError yT^(-2.0)
    @test (1+xT)^(3//2) == ((1+xT)^0.5)^3
    @test real(xH) == xH
    @test imag(xH) == zero(xH)
    @test conj(im*yH) == (im*yH)'
    @test conj(im*yT) == (im*yT)'
    @test real( exp(1im * xT)) == cos(xT)
    @test get_coeff(convert(TaylorN{Rational{Int}},cos(xT)),[4,0]) ==
        1//factorial(4)
    cr = convert(TaylorN{Rational{Int}},cos(xT))
    @test get_coeff(cr,[4,0]) == 1//factorial(4)
    @test imag((exp(yT))^(-1im)') == sin(yT)
    exy = exp( xT+yT )
    @test evaluate(exy) == 1
    @test evaluate(exy,[0.1im,0.01im]) == exp(0.11im)
    @test isapprox(evaluate(exy, [1,1]), e^2)
    txy = tan(xT+yT)
    @test get_coeff(txy,[8,7]) == 929569/99225
    ptxy = xT + yT + (1/3)*( xT^3 + yT^3 ) + xT^2*yT + xT*yT^2
    @test getindex(tan(TaylorN(1)+TaylorN(2)),1:5) == ptxy.coeffs[1:5]
    @test evaluate(xH*yH,[1.0,2.0]) == 2.0
    v = zeros(Int, 2)
    @test evaluate!([xT, yT], ones(Int, 2), v) == nothing
    @test v == ones(2)

    @test evaluate(sin(asin(xT+yT)), [1.0,0.5]) == 1.5
    @test evaluate(asin(sin(xT+yT)), [1.0,0.5]) == 1.5
    @test tan(atan(xT+yT)) == xT+yT
    @test atan(tan(xT+yT)) == xT+yT
    @test asin(xT+yT) + acos(xT+yT) == pi/2

    @test -sinh(xT+yT) + cosh(xT+yT) == exp(-(xT+yT))
    @test  sinh(xT+yT) + cosh(xT+yT) == exp(xT+yT)
    @test evaluate(- sinh(xT+yT)^2 + cosh(xT+yT)^2 , rand(2)) == 1
    @test evaluate(- sinh(xT+yT)^2 + cosh(xT+yT)^2 , zeros(2)) == 1
    @test tanh(xT + yT + 0im) == -1im * tan((xT+yT)*1im)
    @test cosh(xT+yT) == real(cos(im*(xT+yT)))
    @test sinh(xT+yT) == imag(sin(im*(xT+yT)))

    xx = 1.0*zeroT
    TaylorSeries.add!(xx, 1.0*xT, 2yT, 1)
    @test xx[2] == HomogeneousPolynomial([1,2])
    TaylorSeries.subst!(xx, 1.0*xT, yT, 1)
    @test xx[2] == HomogeneousPolynomial([1,-1])
    TaylorSeries.div!(xx, 1.0+xT, 1.0+xT, 0)
    @test xx[1] == 1.0
    TaylorSeries.pow!(xx, 1.0+xT, 1.5, 0)
    @test xx[1] == 1.0
    TaylorSeries.sqrt!(xx, 1.0+xT, 0)
    @test xx[1] == 1.0
    TaylorSeries.exp!(xx, xT, 0)
    @test xx[1] == 1.0
    TaylorSeries.log!(xx, 1.0+xT, 0)
    @test xx[1] == 0.0
    cxx = zero(xx)
    TaylorSeries.sincos!(xx, cxx, 1.0*xT, 0)
    @test xx[1] == 0.0
    @test cxx[1] == 1.0
    TaylorSeries.tan!(xx, 1.0*xT, cxx, 0)
    @test xx[1] == 0.0
    @test cxx[1] == 0.0
    TaylorSeries.asin!(xx, 1.0*xT, cxx, 0)
    @test xx[1] == 0.0
    @test cxx[1] == 1.0
    TaylorSeries.acos!(xx, 1.0*xT, cxx, 0)
    @test xx[1] == acos(0.0)
    @test cxx[1] == 1.0
    TaylorSeries.atan!(xx, 1.0*xT, cxx, 0)
    @test xx[1] == 0.0
    @test cxx[1] == 1.0
    TaylorSeries.sinhcosh!(xx, cxx, 1.0*xT, 0)
    @test xx[1] == 0.0
    @test cxx[1] == 1.0
    TaylorSeries.tanh!(xx, 1.0*xT, cxx, 0)
    @test xx[1] == 0.0
    @test cxx[1] == 0.0

    g1(xT,yT) = xT^3 + 3yT^2 - 2xT^2 * yT - 7xT + 2
    g2(xT,yT) = yT + xT^2 - xT^4
    f1 = g1(xT,yT)
    f2 = g2(xT,yT)
    @test gradient(f1) == [ 3*xT^2-4*xT*yT-TaylorN(7,0), 6*yT-2*xT^2 ]
    @test ∇(f2) == [2*xT - 4*xT^3, TaylorN(1,0)]
    @test jacobian([f1,f2], [2,1]) == jacobian( [g1(xT+2,yT+1), g2(xT+2,yT+1)] )
    jac = Array{Int64}(2, 2)
    jacobian!(jac, [g1(xT+2,yT+1), g2(xT+2,yT+1)])
    @test jac == jacobian( [g1(xT+2,yT+1), g2(xT+2,yT+1)] )
    jacobian!(jac, [f1,f2], [2,1])
    @test jac == jacobian([f1,f2], [2,1])
    @test [xT yT]*hessian(f1*f2)*[xT, yT] == [ 2*TaylorN((f1*f2)[3]) ]
    @test hessian(f1^2)/2 == [ [49,0] [0,12] ]
    @test hessian(f1-f2-2*f1*f2) == (hessian(f1-f2-2*f1*f2))'
    @test hessian(f1-f2,[1,-1]) == hessian(g1(xT+1,yT-1)-g2(xT+1,yT-1))
    hes = Array{Int64}(2, 2)
    hessian!(hes, f1*f2)
    @test hes == hessian(f1*f2)
    @test [xT yT]*hes*[xT, yT] == [ 2*TaylorN((f1*f2)[3]) ]
    hessian!(hes, f1^2)
    @test hes/2 == [ [49,0] [0,12] ]
    hessian!(hes, f1-f2-2*f1*f2)
    @test hes == hes'
    hes1 = hes2 = Array{Int64}(2, 2)
    hessian!(hes1,f1-f2,[1,-1])
    hessian!(hes2,g1(xT+1,yT-1)-g2(xT+1,yT-1))
    @test hes1 == hes2

    @test string(-xH) == " - 1 x₁"
    @test string(xT^2) == " 1 x₁² + 𝒪(‖x‖¹⁸)"
    @test string(1im*yT) == " ( 1 im ) x₂ + 𝒪(‖x‖¹⁸)"
    @test string(xT-im*yT) == "  ( 1 ) x₁ - ( 1 im ) x₂ + 𝒪(‖x‖¹⁸)"

    @test_throws ArgumentError abs(xT)
    @test_throws AssertionError 1/x
    @test_throws AssertionError zero(x)/zero(x)
    @test_throws ArgumentError sqrt(x)
    @test_throws AssertionError x^(-2)
    @test_throws ArgumentError log(x)
    @test_throws AssertionError cos(x)/sin(y)
    @test_throws BoundsError xH[20]
    @test_throws BoundsError xT[20]
end
