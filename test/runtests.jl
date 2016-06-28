# This file is part of TaylorSeries.jl, MIT licensed
#
# Tests for TaylorSeries implementation
using TaylorSeries
using FactCheck
using Compat
import Compat.String

FactCheck.setstyle(:compact)
# FactCheck.onlystats(true)

facts("Tests for Taylor1 expansions") do
    ta(a) = Taylor1([a,one(a)],15)
    t = Taylor1(Int,15)
    tim = im*t
    zt = zero(t)
    ot = 1.0*one(t)
    tol1 = eps(1.0)

    @fact Taylor1([0,1,0,0]) == Taylor1(3)  --> true
    @fact get_coeff(Taylor1(Complex128,3),1) == complex(1.0,0.0) --> true
    @fact eltype(convert(Taylor1{Complex128},ot)) == Complex128  --> true
    @fact eltype(convert(Taylor1{Complex128},1)) == Complex128  --> true
    @fact convert(Taylor1{Complex{Int}},[0,2]) == (2+0im)*t  --> true
    @fact convert(Taylor1{BigFloat},[0.0, 1.0]) == ta(big(0.0))  --> true
    @fact promote(0,Taylor1(1.0,0)) == (zt,ot)  --> true
    @fact eltype(promote(ta(0.0),zeros(Int,2))[2]) == Float64  --> true
    @fact eltype(promote(0,Taylor1(ot))[1]) == Float64  --> true
    @fact eltype(promote(1.0+im, zt)[1]) == Complex{Float64}  --> true
    @fact eltype(TaylorSeries.fixshape(zt,ot)[1]) == Float64  --> true

    @fact length(Taylor1(10)) == 10  --> true
    @fact length(TaylorSeries.fixshape(zt,Taylor1([1.0]))[2]) == 15  --> true
    @fact eltype(TaylorSeries.fixshape(zt,Taylor1([1.0]))[1]) == Float64  --> true
    @fact TaylorSeries.firstnonzero(t) == 1  --> true
    @fact TaylorSeries.firstnonzero(zt) == zt.order+1  --> true

    @fact t == Taylor1(ta(0),15)  --> true
    @fact ot == 1  --> true
    @fact 0.0 == zt  --> true
    @fact get_coeff(tim,1) == complex(0,1)  --> true
    @fact zt+1.0 == ot  --> true
    @fact 1.0-ot == zt  --> true
    @fact t+t == 2t  --> true
    @fact t-t == zt  --> true
    @fact +t == -(-t)  --> true

    tsquare = Taylor1([0,0,1],15)
    @fact t * true == t  --> true
    @fact false * t == zero(t)  --> true
    @fact t^0 == t^0.0 == one(t)  --> true
    @fact t*t == tsquare  --> true
    @fact t*1 == t -->  true
    @fact 0*t == zt -->  true
    @fact (-t)^2 == tsquare  --> true
    @fact t^3 == tsquare*t  --> true
    @fact tsquare/t == t  --> true
    @fact t/(t*3) == (1/3)*ot  --> true
    @fact t/3im == -tim/3  --> true
    @fact 1/(1-t) == Taylor1(ones(t.order+1))  --> true
    @fact Taylor1([0,1,1])/t == t+1  --> true
    @fact (t+im)^2 == tsquare+2im*t-1  --> true
    @fact (t+im)^3 == Taylor1([-1im,-3,3im,1],15)  --> true
    @fact (t+im)^4 == Taylor1([1,-4im,-6,4im,1],15)  --> true
    @fact imag(tsquare+2im*t-1) == 2t  --> true
    @fact (Rational(1,2)*tsquare).coeffs[3] == 1//2  --> true
    @fact t^2/tsquare == ot  --> true
    @fact ((1+t)^(1/3)).coeffs[3]+1/9 <= tol1  --> true
    @fact 1-tsquare == (1+t)-t*(1+t)  --> true
    @fact (1-tsquare)^2 == (1+t)^2.0 * (1-t)^2.0  --> true
    @fact (sqrt(1+t)).coeffs[3] == -1/8  --> true
    @fact ((1-tsquare)^(1//2))^2 == 1-tsquare  --> true
    @fact ((1-t)^(1//4)).coeffs[15] == -4188908511//549755813888  --> true
    @fact abs(((1+t)^3.2).coeffs[14] + 5.4021062656e-5) < tol1  --> true

    @fact isapprox( rem(4.1 + t,4).coeffs[1], (0.1 + t).coeffs[1] )  --> true
    @fact isapprox( mod(4.1 + t,4).coeffs[1], (0.1 + t).coeffs[1] )  --> true
    @fact isapprox( mod2pi(2pi+0.1+t).coeffs[1],(0.1 + t).coeffs[1])  --> true

    @fact abs(ta(1))  --> ta(1)
    @fact abs(ta(-1.0))  --> -ta(-1.0)

    @fact log(exp(tsquare)) == tsquare  --> true
    @fact exp(log(1-tsquare)) == 1-tsquare  --> true
    @fact log((1-t)^2) == 2*log(1-t)  --> true
    @fact real(exp(tim)) == cos(t)  --> true
    @fact imag(exp(tim)) == sin(t)  --> true
    @fact exp(conj(tim)) == cos(t)-im*sin(t) == exp(tim')  --> true
    @fact (exp(t))^(2im) == cos(2t)+im*sin(2t)  --> true
    @fact (exp(t))^Taylor1([-5.2im]) == cos(5.2t)-im*sin(5.2t)  --> true
    @fact get_coeff(convert(Taylor1{Rational{Int}},cos(t)),8) ==
        1//factorial(8)  --> true
    @fact abs((tan(t)).coeffs[8]- 17/315) < tol1  --> true
    @fact abs((tan(t)).coeffs[14]- 21844/6081075) < tol1  --> true
    @fact evaluate(exp(Taylor1([0,1],17)),1.0) == 1.0*e  --> true
    @fact evaluate(exp(Taylor1([0,1],1))) == 1.0  --> true
    @fact evaluate(exp(t),t^2) == exp(t^2)  --> true

    @fact derivative(5, exp(ta(1.0))) == exp(1.0)  --> true
    @fact derivative(3, exp(ta(1.0pi))) == exp(1.0pi)  --> true
    @fact isapprox(derivative(10, exp(ta(1.0pi))) , exp(1.0pi) )  --> true
    @fact integrate(derivative(exp(t)),1) == exp(t)  --> true
    @fact integrate(cos(t)) == sin(t)  --> true

    @fact promote(ta(0.0), t) == (ta(0.0),ta(0.0))  --> true

    @fact_throws ArgumentError abs(ta(big(0)))
    @fact_throws ArgumentError 1/t
    @fact_throws ArgumentError zt/zt
    @fact_throws ArgumentError t^1.5
    @fact_throws DomainError t^(-2)
    @fact_throws ArgumentError sqrt(t)
    @fact_throws ArgumentError log(t)
    @fact_throws ArgumentError cos(t)/sin(t)
    @fact_throws AssertionError derivative(30, exp(ta(1.0pi)))

    @fact string(ta(-3)) == " - 3 + 1 t + 𝒪(t¹⁶)"  --> true
    @fact TaylorSeries.pretty_print(ta(3im)) ==
        " ( 3 im )  + ( 1 ) t + 𝒪(t¹⁶)"  --> true
end

facts("Tests for HomogeneousPolynomial and TaylorN") do

    @fact eltype(set_variables(Int, "x", numvars=2, order=6))  --> TaylorN{Int}
    @fact eltype(set_variables("x", numvars=2, order=6))  --> TaylorN{Float64}
    @fact eltype(set_variables(BigInt, "x y", order=6))  --> TaylorN{BigInt}
    @fact eltype(set_variables("x y", order=6))  --> TaylorN{Float64}

    @fact TaylorSeries.coeff_table[2][1] == [1,0]  --> true
    @fact TaylorSeries.index_table[2][1] == 7  --> true
    @fact TaylorSeries.in_base(get_order(),[2,1]) == 15  --> true
    @fact TaylorSeries.pos_table[4][15] == 2  --> true

    @fact get_order() == 6  --> true
    @fact get_numvars() == 2  --> true

    x, y = set_variables("x y", order=6)
    @fact x.order == 6 --> true
    @fact TaylorSeries.set_variable_names(["x","y"]) == ["x", "y"]  --> true
    @fact TaylorSeries.get_variable_names() == ["x", "y"]  --> true

    set_variables("x", numvars=2, order=17)

    xH = HomogeneousPolynomial([1,0])
    yH = HomogeneousPolynomial([0,1],1)
    xT = TaylorN(xH,17)
    yT = taylorN_variable(Int64, 2, 17)
    zeroT = zero( TaylorN([xH],1) )
    uT = one(convert(TaylorN{Float64},yT))
    @fact ones(xH,1) == [HomogeneousPolynomial(1), xH+yH]  --> true
    @fact ones(HomogeneousPolynomial{Complex{Int}},0) ==
        [HomogeneousPolynomial(1+0im)]  --> true
    @fact get_order(zeroT) == 1  --> true
    @fact get_coeff(xT,[1,0]) == 1  --> true
    @fact get_coeff(yH,[1,0]) == 0  --> true
    @fact convert(HomogeneousPolynomial{Int64},[1,1]) == xH+yH  --> true
    @fact convert(HomogeneousPolynomial{Float64},[2,-1]) == 2.0xH-yH  --> true
    @fact convert(TaylorN{Float64}, yH) == 1.0*yT  --> true
    @fact convert(TaylorN{Float64}, [xH,yH]) == xT+1.0*yT  --> true
    @fact convert(TaylorN{Int}, [xH,yH]) == xT+yT  --> true
    @fact promote(xH, [1,1])[2] == xH+yH  --> true
    @fact promote(xH, yT)[1] == xT  --> true
    @fact promote(xT, [xH,yH])[2] == xT+yT  --> true
    @fact typeof(promote(im*xT,[xH,yH])[2]) == TaylorN{Complex{Int64}}  --> true
    @fact TaylorSeries.fixorder(taylorN_variable(1,1),17) == xT  --> true
    @fact TaylorSeries.iszero(zeroT.coeffs[2])  -->  true

    @fact HomogeneousPolynomial(xH,1) == HomogeneousPolynomial(xH)  --> true
    @fact eltype(xH) == Int  --> true
    @fact length(xH) == 2  --> true
    @fact zero(xH) == 0*xH  --> true
    @fact one(yH) == xH+yH  --> true
    @fact xH * true == xH  --> true
    @fact false * yH == zero(yH)  --> true
    @fact get_order(yH) == 1  --> true
    @fact get_order(xT) == 17  --> true
    @fact xT * true == xT  --> true
    @fact false * yT == zero(yT)  --> true

    @fact xT == TaylorN([xH])  --> true
    @fact one(xT) == TaylorN(1,5)  --> true
    @fact TaylorN(zeroT,5) == 0  --> true
    @fact TaylorN(uT) == convert(TaylorN{Complex},1)  --> true
    @fact get_numvars() == 2  --> true
    @fact length(uT) == get_order()+1  --> true
    @fact eltype(convert(TaylorN{Complex128},1)) == Complex128  --> true

    @fact 1+xT+yT == TaylorN(1,1) + TaylorN([xH,yH],1)  --> true
    @fact xT-yT-1 == TaylorN([-1,xH-yH])  --> true
    @fact xT*yT == TaylorN([HomogeneousPolynomial([0,1,0],2)])  --> true
    @fact (1/(1-xT)).coeffs[4] == HomogeneousPolynomial(1.0,3)  --> true
    @fact xH^20 == HomogeneousPolynomial([0],get_order())  --> true
    @fact (yT/(1-xT)).coeffs[5] == xH^3 * yH  --> true
    @fact mod(1+xT,1) == +xT  --> true
    @fact (rem(1+xT,1)).coeffs[1] == 0  --> true
    @fact abs(1-xT)  --> 1-xT
    @fact abs(-1-xT)  --> 1+xT
    @fact derivative(mod2pi(2pi+yT^3),2) == derivative(yT^3,2)  --> true
    @fact derivative(yT) == zeroT  --> true
    @fact -xT/3im == im*xT/3  --> true
    @fact (xH/3im)' == im*xH/3  --> true

    @fact derivative(2xT*yT^2,1) == 2yT^2  --> true
    @fact xT*xT^3 == xT^4  --> true
    txy = 1.0 + xT*yT - 0.5*xT^2*yT + (1/3)*xT^3*yT + 0.5*xT^2*yT^2
    @fact (1+taylorN_variable(1,4))^taylorN_variable(2,4) == txy  --> true
    @fact_throws DomainError yT^(-2)
    @fact_throws DomainError yT^(-2.0)
    @fact (1+xT)^(3//2) == ((1+xT)^0.5)^3  --> true
    @fact real(xH) == xH  --> true
    @fact imag(xH) == zero(xH)  --> true
    @fact conj(im*yH) == (im*yH)'  --> true
    @fact conj(im*yT) == (im*yT)'  --> true
    @fact real( exp(1im * xT)) == cos(xT)  --> true
    @fact get_coeff(convert(TaylorN{Rational{Int}},cos(xT)),[4,0]) ==
        1//factorial(4)  --> true
    cr = convert(TaylorN{Rational{Int}},cos(xT))
    @fact get_coeff(cr,[4,0]) == 1//factorial(4)  --> true
    @fact imag((exp(yT))^(-1im)') == sin(yT)  --> true
    exy = exp( xT+yT )
    @fact evaluate(exy) == 1  --> true
    @fact evaluate(exy,[0.1im,0.01im]) == exp(0.11im)  --> true
    @fact isapprox(evaluate(exy, [1,1]), e^2)  --> true
    txy = tan(xT+yT)
    @fact get_coeff(txy,[8,7]) == 929569/99225  --> true
    ptxy = xT + yT + (1/3)*( xT^3 + yT^3 ) + xT^2*yT + xT*yT^2
    @fact tan(taylorN_variable(1,4)+taylorN_variable(2,4)) == ptxy  --> true
    @fact evaluate(xH*yH,[1.0,2.0]) == 2.0  --> true

    g1(xT,yT) = xT^3 + 3yT^2 - 2xT^2 * yT - 7xT + 2
    g2(xT,yT) = yT + xT^2 - xT^4
    f1 = g1(xT,yT)
    f2 = g2(xT,yT)
    @fact gradient(f1) == [ 3*xT^2-4*xT*yT-TaylorN(7), 6*yT-2*xT^2 ]  --> true
    @fact ∇(f2) == [2*xT - 4*xT^3, TaylorN(1)]  --> true
    @fact jacobian([f1,f2], [2,1]) ==
        jacobian( [g1(xT+2,yT+1), g2(xT+2,yT+1)] )  --> true
    @fact [xT yT]*hessian(f1*f2)*[xT, yT] ==
        [ 2*TaylorN((f1*f2).coeffs[3]) ]  --> true
    @fact hessian(f1^2)/2 == [ [49,0] [0,12] ]  --> true
    @fact hessian(f1-f2-2*f1*f2) == (hessian(f1-f2-2*f1*f2))'  --> true
    @fact hessian(f1-f2,[1,-1]) == hessian(g1(xT+1,yT-1)-g2(xT+1,yT-1))  --> true

    @fact string(-xH) == " - 1 x₁"  --> true
    @fact string(xT^2) == " 1 x₁² + 𝒪(‖x‖¹⁸)"  --> true
    @fact string(1im*yT) == " ( 1 im ) x₂ + 𝒪(‖x‖¹⁸)"  --> true
    @fact string(xT-im*yT) == "  ( 1 ) x₁ - ( 1 im ) x₂ + 𝒪(‖x‖¹⁸)"  --> true

    @fact_throws ArgumentError abs(xT)
    @fact_throws AssertionError 1/x
    @fact_throws AssertionError zero(x)/zero(x)
    @fact_throws ArgumentError sqrt(x)
    @fact_throws AssertionError x^(-2)
    @fact_throws ArgumentError log(x)
    @fact_throws AssertionError cos(x)/sin(y)
end

facts("Testing an identity proved by Euler (8 variables)") do
    make_variable(name, index::Int) = string(name, TaylorSeries.subscriptify(index))

    @compat variable_names = String[make_variable("α", i) for i in 1:4]
    append!(variable_names, [make_variable("β", i) for i in 1:4])

    a1, a2, a3, a4, b1, b2, b3, b4 = set_variables(variable_names, order=4)

    lhs1 = a1^2 + a2^2 + a3^2 + a4^2
    lhs2 = b1^2 + b2^2 + b3^2 + b4^2
    lhs = lhs1 * lhs2

    rhs1 = (a1*b1 - a2*b2 - a3*b3 - a4*b4)^2
    rhs2 = (a1*b2 + a2*b1 + a3*b4 - a4*b3)^2
    rhs3 = (a1*b3 - a2*b4 + a3*b1 + a4*b2)^2
    rhs4 = (a1*b4 + a2*b3 - a3*b2 + a4*b1)^2

    rhs = rhs1 + rhs2 + rhs3 + rhs4

    @fact lhs == rhs  --> true
    v = randn(8)
    @fact evaluate( rhs, v) == evaluate( lhs, v)  --> true
end

facts("High order polynomials test inspired by Fateman (takes a few seconds))") do
    x, y, z, w = set_variables(Int128, "x", numvars=4, order=40)

    function fateman2(degree::Int)
        T = Int128
        oneH = HomogeneousPolynomial(one(T), 0)
        # s = 1 + x + y + z + w
        s = TaylorN(
         [oneH, HomogeneousPolynomial([one(T),one(T),one(T),one(T)],1)], degree)
        s = s^degree
        # s is converted to order 2*ndeg
        s = TaylorN(s, 2*degree)
        return s^2 + s
    end

    function fateman3(degree::Int)
        s = x + y + z + w + 1
        s = s^degree
        s * (s+1)
    end

    f2 = fateman2(20)
    f3 = fateman3(20)
    c = get_coeff(f2,[1,6,7,20])
    @fact c == 128358585324486316800 --> true
    @fact get_coeff(f2,[1,6,7,20]) == c --> true
end

facts("Matrix multiplication for Taylor1") do
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
        @fact A*B==Y -->true
        # Y should be extended after the multilpication
        @fact reduce(&, [y1.order for y1 in Y] .== Y[1].order) --> true
        # B should be unchanged
        @fact B==Bcopy --> true

        # is the result compatible with the matrix multiplication?  We
        # only check the zeroth order of the Taylor series.
        y1=sum(Y).coeffs[1]
        Y=A*B1[:,1]
        y2=sum(Y)

        # There is a small numerical error when comparing the generic
        # multiplication and the specialized version
        @fact abs(y1-y2) < n1*(eps(y1)+eps(y2)) --> true

        @fact_throws DimensionMismatch A_mul_B!(Y,A[:,1:end-1],B)
        @fact_throws DimensionMismatch A_mul_B!(Y,A[1:end-1,:],B)
        @fact_throws DimensionMismatch A_mul_B!(Y,A,B[1:end-1])
        @fact_throws DimensionMismatch A_mul_B!(Y[1:end-1],A,B)
    end

end


exitstatus()
