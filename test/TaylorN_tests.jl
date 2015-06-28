using TaylorSeries
using FactCheck

facts("Tests for HomogeneousPolynomial and TaylorN") do

    x, y = set_variables("x", numvars=2, order=6)

    @fact get_order() == 6  => true
    @fact get_numvars() == 2  => true

    @fact x.order == 6 => true


    set_variables("x", numvars=2, order=17)

    xH = HomogeneousPolynomial([1,0])
    yH = HomogeneousPolynomial([0,1],1)

    xT = TaylorN(xH,17)
    yT = taylorN_variable(Int64, 2, 17)

    zeroT = zero( TaylorN([xH],1) )
    uT = one(convert(TaylorN{Float64},yT))

    @fact ones(xH,1) == [HomogeneousPolynomial(1), xH+yH]  => true
    @fact ones(HomogeneousPolynomial{Complex{Int}},0) ==
        [HomogeneousPolynomial(1+0im)]  => true
    @fact zeroT.order == 1  => true
    @fact get_coeff(xT, (1,0) ) == 1  => true
    @fact get_coeff(yH, (1,0) ) == 0  => true
    @fact convert(HomogeneousPolynomial{Int64}, [1,1]) == xH+yH  => true
    @fact convert(HomogeneousPolynomial{Float64}, [2,-1]) == 2.0xH-yH  => true
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
    @fact get_order(yH) == 1  => true
    @fact get_order(xT) == 17  => true

    @fact TaylorN(zeroT,5) == 0  => true
    @fact TaylorN(uT) == convert(TaylorN{Complex},1)  => true
    @fact get_numvars() == 2  => true
    @fact length(uT) == get_order()+1  => true
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
    exy = exp( xT+yT )
    @fact evaluate(exy) == 1  => true
    @fact isapprox(evaluate(exy, [1,1]), e^2)  => true
    @fact evaluate(exy) == 1  => true
    txy = tan(xT+yT)
    @fact get_coeff(txy, (8,7) ) == 929569/99225  => true
    ptxy = xT + yT + (1/3)*( xT^3 + yT^3 ) + xT^2*yT + xT*yT^2
    @fact tan(taylorN_variable(1,4)+taylorN_variable(2,4)) == ptxy  => true

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
