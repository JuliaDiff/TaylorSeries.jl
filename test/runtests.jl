# Tests for TaylorSeries implementation
using TaylorSeries
using FactCheck
using Compat

FactCheck.setstyle(:compact)
# FactCheck.onlystats(true)

facts("Tests for Taylor1 expansions") do
    ta(a) = Taylor1([a,one(a)],15)
    t = taylor1_variable(Int,15)
    tim = im*t
    zt = zero(t)
    ot = 1.0*one(t)
    tol1 = eps(1.0)

    @fact Taylor1([0,1,0,0]) == taylor1_variable(3)  => true
    @fact eltype(convert(Taylor1{Complex128},ot)) == Complex128  => true
    @fact eltype(convert(Taylor1{Complex128},1)) == Complex128  => true
    @fact convert(Taylor1{Complex{Int}},[0,2]) == (2+0im)*t  => true
    @fact promote(0,Taylor1(1.0,0)) == (zt,ot)  => true
    @fact eltype(promote(0,Taylor1(ot))[1]) == Float64  => true
    @fact eltype(promote(1.0+im, zt)[1]) == Complex{Float64}  => true
    @fact eltype(TaylorSeries.fixshape(zt,ot)[1]) == Float64  => true

    @fact length(Taylor1(0)) == 0  => true
    @fact length(TaylorSeries.fixshape(zt,convert(Taylor1{Int64},[0]))[1]) ==
        15  => true
    @fact TaylorSeries.firstnonzero(t) == 1  => true
    @fact TaylorSeries.firstnonzero(zt) == zt.order+1  => true

    @fact t == Taylor1(ta(0),15)  => true
    @fact ot == 1  => true
    @fact 0.0 == zt  => true
    @fact get_coeff(t,1) == 1  => true
    @fact zt+1 == ot  => true
    @fact t+t == 2t  => true
    @fact t-t == zt  => true

    tsquare = Taylor1([0,0,1],15)
    @fact t^0 == t^0.0 == one(t)  => true
    @fact t*t == tsquare  => true
    @fact (-t)^2 == tsquare  => true
    @fact t^3 == tsquare*t  => true
    @fact tsquare/t == t  => true
    @fact t/(t*3) == (1/3)*ot  => true
    @fact t/3im == -tim/3  => true
    @fact 1/(1-t) == Taylor1(ones(t.order+1))  => true
    @fact Taylor1([0,1,1])/t == t+1  => true
    @fact (t+im)^2 == tsquare+2im*t-1  => true
    @fact (t+im)^3 == Taylor1([-1im,-3,3im,1],15)  => true
    @fact (t+im)^4 == Taylor1([1,-4im,-6,4im,1],15)  => true
    @fact imag(tsquare+2im*t-1) == 2t  => true
    @fact (Rational(1,2)*tsquare).coeffs[3] == 1//2  => true
    @fact t^2/tsquare == ot  => true
    @fact ((1+t)^(1/3)).coeffs[3]+1/9 <= tol1  => true
    @fact 1-tsquare == (1+t)-t*(1+t)  => true
    @fact (1-tsquare)^2 == (1+t)^2.0 * (1-t)^2.0  => true
    @fact (sqrt(1+t)).coeffs[3] == -1/8  => true
    @fact ((1-tsquare)^(1//2))^2 == 1-tsquare  => true
    @fact ((1-t)^(1//4)).coeffs[15] == -4188908511//549755813888  => true
    @fact abs(((1+t)^3.2).coeffs[14] + 5.4021062656e-5) < tol1  => true

    @fact isapprox( rem(4.1 + t,4).coeffs[1], (0.1 + t).coeffs[1] )  => true
    @fact isapprox( mod(4.1 + t,4).coeffs[1], (0.1 + t).coeffs[1] )  => true
    @fact isapprox( mod2pi(2pi+0.1+t).coeffs[1],(0.1 + t).coeffs[1])  => true

    @fact log(exp(tsquare)) == tsquare  => true
    @fact exp(log(1-tsquare)) == 1-tsquare  => true
    @fact log((1-t)^2) == 2*log(1-t)  => true
    @fact real(exp(tim)) == cos(t)  => true
    @fact imag(exp(tim)) == sin(t)  => true
    @fact exp(tim') == cos(t)-im*sin(t)  => true
    @fact (exp(t))^(2im) == cos(2t)+im*sin(2t)  => true
    @fact (exp(t))^Taylor1(-5.2im) == cos(5.2t)-im*sin(5.2t)  => true
    @fact abs((tan(t)).coeffs[8]- 17/315) < tol1  => true
    @fact abs((tan(t)).coeffs[14]- 21844/6081075) < tol1  => true
    @fact evalTaylor(exp(Taylor1([0,1],17)),1.0) == 1.0*e  => true
    @fact evalTaylor(exp(Taylor1([0,1],1))) == 1.0  => true
    @fact evalTaylor(exp(t),t^2) == exp(t^2)  => true

    @fact deriv( exp(ta(1.0)), 5 ) == exp(1.0)  => true
    @fact deriv( exp(ta(1.0pi)), 3 ) == exp(1.0pi)  => true
    @fact isapprox( deriv(exp(ta(1.0pi)), 10) , exp(1.0pi) )  => true
    @fact integTaylor(diffTaylor(exp(t)),1) == exp(t)  => true
    @fact integTaylor(cos(t)) == sin(t)  => true

    @fact promote(ta(0.0), t) == (ta(0.0),ta(0.0))  => true

    @fact_throws ErrorException 1/t
    @fact_throws ErrorException zt/zt
    @fact_throws ErrorException t^1.5
    @fact_throws DomainError t^(-2)
    @fact_throws ErrorException sqrt(t)
    @fact_throws ErrorException log(t)
    @fact_throws ErrorException cos(t)/sin(t)
    # @fact_throws AssertionError deriv( exp(ta(1.0pi)), 30 )
end

facts("Tests for HomogeneousPolynomial and TaylorN") do
    set_params_TaylorN(6,2)
    @fact set_numVars(2) == (6,2)  => true
    @fact set_maxOrder(6) == (6,2)  => true
    @fact get_maxOrder() == 6  => true
    @fact get_numVars() == 2  => true

    set_params_TaylorN(17,2)
    xH = HomogeneousPolynomial([1,0])
    yH = HomogeneousPolynomial([0,1],1)
    xT = TaylorN(xH,17)
    yT = taylorN_variable(Int64, 2, 17)
    zeroT = zero( TaylorN([xH],1) )
    uT = one(convert(TaylorN{Float64},yT))
    @fact ones(xH,1) == [HomogeneousPolynomial(1), xH+yH]  => true
    @fact ones(HomogeneousPolynomial{Complex{Int}},0) ==
        [HomogeneousPolynomial(1+0im)]  => true
    @fact get_maxOrder(zeroT) == 1  => true
    @fact get_coeff(xT,[1,0]) == 1  => true
    @fact get_coeff(yH,[1,0]) == 0  => true
    @fact convert(HomogeneousPolynomial{Int64},[1,1]) == xH+yH  => true
    @fact convert(HomogeneousPolynomial{Float64},[2,-1]) == 2.0xH-yH  => true
    @fact convert(TaylorN{Float64}, yH) == 1.0*yT  => true
    @fact convert(TaylorN{Float64}, [xH,yH]) == xT+1.0*yT  => true
    @fact convert(TaylorN{Int}, [xH,yH]) == xT+yT  => true
    @fact TaylorSeries.fixorder(taylorN_variable(1,1),17) == xT  => true
    @fact TaylorSeries.iszero(zeroT.coeffs[2])  =>  true

    @fact HomogeneousPolynomial(xH,1) == HomogeneousPolynomial(xH)  => true
    @fact eltype(xH) == Int  => true
    @fact length(xH) == 2  => true
    @fact zero(xH) == 0*xH  => true
    @fact one(yH) == xH+yH  => true
    @fact get_maxOrder(yH) == 1  => true
    @fact get_maxOrder(xT) == 17  => true

    @fact TaylorN(zeroT,5) == 0  => true
    @fact TaylorN(uT) == convert(TaylorN{Complex},1)  => true
    @fact get_numVars() == 2  => true
    @fact length(uT) == get_maxOrder()+1  => true
    @fact eltype(convert(TaylorN{Complex128},1)) == Complex128  => true

    @fact 1+xT+yT == TaylorN(1,1) + TaylorN([xH,yH],1)  => true
    @fact xT-yT == TaylorN([xH-yH])  => true
    @fact xT*yT == TaylorN([HomogeneousPolynomial([0,1,0],2)])  => true
    @fact (1/(1-xT)).coeffs[4] == HomogeneousPolynomial(1.0,3)  => true
    @fact (yT/(1-xT)).coeffs[5] == xH^3 * yH  => true
    @fact mod(1+xT,1) == +xT  => true
    @fact (rem(1+xT,1)).coeffs[1] == 0  => true
    @fact diffTaylor(mod2pi(2pi+yT^3),2) == diffTaylor(yT^3,2)  => true
    @fact diffTaylor(yT) == zeroT  => true
    @fact -xT/3im == im*xT/3  => true
    @fact (xH/3im)' == im*xH/3  => true

    @fact diffTaylor(2xT*yT^2,1) == 2yT^2  => true
    @fact xT*xT^3 == xT^4  => true
    txy = 1.0 + xT*yT - 0.5*xT^2*yT + (1/3)*xT^3*yT + 0.5*xT^2*yT^2
    @fact (1+taylorN_variable(1,4))^taylorN_variable(2,4) == txy  => true
    @fact_throws DomainError yT^(-2)
    @fact_throws DomainError yT^(-2.0)
    @fact (1+xT)^(3//2) == ((1+xT)^0.5)^3  => true
    @fact real( exp(1im * xT)) == cos(xT)  => true
    @fact imag((exp(yT))^(-1im)') == sin(yT)  => true
    @fact evalTaylor(exp( xT+yT )) == 1  => true
    @fact isapprox(evalTaylor(exp( xT+yT ), [1,1]), e^2)  => true
    txy = xT + yT + (1/3)*( xT^3 + yT^3 ) + xT^2*yT + xT*yT^2
    @fact tan(taylorN_variable(1,4)+taylorN_variable(2,4)) == txy  => true

    g1(xT,yT) = xT^3 + 3yT^2 - 2xT^2 * yT - 7xT + 2
    g2(xT,yT) = yT + xT^2 - xT^4
    f1 = g1(xT,yT)
    f2 = g2(xT,yT)
    @fact gradient(f1) == [ 3*xT^2-4*xT*yT-TaylorN(7), 6*yT-2*xT^2 ]  => true
    @fact ∇(f2) == [2*xT - 4*xT^3, TaylorN(1)]  => true
    @fact jacobian([f1,f2], [2,1]) ==
        jacobian( [g1(xT+2,yT+1), g2(xT+2,yT+1)] )  => true
    @fact [xT yT]*hessian(f1*f2)*[xT, yT] ==
        [ 2*TaylorN((f1*f2).coeffs[3]) ]  => true
    @fact hessian(f1^2)/2 == [ [49,0] [0,12] ]  => true
    @fact hessian(f1-f2-2*f1*f2) == (hessian(f1-f2-2*f1*f2))'  => true
    @fact hessian(f1-f2,[1,-1]) == hessian(g1(xT+1,yT-1)-g2(xT+1,yT-1))  => true
end

facts("Testing an identity proved by Euler (8 variables)") do
    @fact set_params_TaylorN(4,8) == (4,8)  => true  # order 4, 8 variables
    # This creates symbols :a1, :b1, ... :a4, :b4 which are independent variables
    for i=1:4
        ai = symbol(string("a",i))
        bi = symbol(string("b",i))
        @eval ($ai) = taylorN_variable(Int,$i,4)
        @eval ($bi) = taylorN_variable(Int,4+($i),4)
    end
    expr_lhs1 = a1^2 + a2^2 + a3^2 + a4^2
    expr_lhs2 = b1^2 + b2^2 + b3^2 + b4^2
    lhs = expr_lhs1 * expr_lhs2
    expr_rhs1 = (a1*b1 - a2*b2 - a3*b3 - a4*b4)^2
    expr_rhs2 = (a1*b2 + a2*b1 + a3*b4 - a4*b3)^2
    expr_rhs3 = (a1*b3 - a2*b4 + a3*b1 + a4*b2)^2
    expr_rhs4 = (a1*b4 + a2*b3 - a3*b2 + a4*b1)^2
    rhs = expr_rhs1 + expr_rhs2 + expr_rhs3 + expr_rhs4
    @fact lhs == rhs  => true
end

exitstatus()
