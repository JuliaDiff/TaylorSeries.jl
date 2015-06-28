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
    @fact evaluate(exp(Taylor1([0,1],17)),1.0) == 1.0*e  => true
    @fact evaluate(exp(Taylor1([0,1],1))) == 1.0  => true
    @fact evaluate(exp(t),t^2) == exp(t^2)  => true

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
