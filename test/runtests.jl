# Simple tests for TaylorSeries implementation
using TaylorSeries
using Base.Test

# Tests for 1-d Tylor expansions
x0 = Taylor(0)
@test length(x0) == 0

x(a) = Taylor([a,one(a)],15)
x0 = x(0)
xI = im*x0
z = zero(x0)
u = 1.0*one(x0)
tol1 = eps(1.0)

@test eltype(convert(Taylor{Complex128},u)) == Complex128
@test eltype(promote(z,Taylor(u))[1]) == Float64
@test eltype(promote(1.0+im, z)[1]) == Complex{Float64}
@test eltype(TaylorSeries.fixshape(z,u)[1]) == Float64
@test length(TaylorSeries.fixshape(Taylor(0,5),z)[1]) == 15
@test TaylorSeries.firstnonzero(x0) == 1
@test TaylorSeries.firstnonzero(z) == z.order+1

@test u == 1
@test 0.0 == z
@test x0.coeffs[2] == 1
@test z+1 == u
@test x0+x0 == 2x0
@test x0-x0 == z

xsquare = Taylor([0,0,1],15)
@test x0*x0 == xsquare
@test (-x0)^2 == xsquare
@test xsquare/x0 == x0
@test x0/(x0*3) == (1/3)*u
@test x0/3im == -xI/3
@test 1/(1-x0) == Taylor(ones(x0.order+1))
@test Taylor([0,1,1])/x0 == x0+1
@test (x0+im)^2 == xsquare+2im*x0-1
@test imag(xsquare+2im*x0-1) == 2x0
@test (Rational(1,2)*xsquare).coeffs[3] == 1//2
@test x0^2/xsquare == u
@test ((1+x0)^(1/3)).coeffs[3]+1/9 <= tol1
@test 1-xsquare == (1+x0)-x0*(1+x0)
@test (1-xsquare)^2 == (1+x0)^2 * (1-x0)^2
@test (sqrt(1+x0)).coeffs[3] == -1/8
@test ((1-x0)^(1/4)).coeffs[15] == -4188908511/549755813888
@test abs(((1+x0)^3.2).coeffs[14] + 5.4021062656e-5) < tol1

@test log(exp(xsquare)) == xsquare
@test exp(log(1-xsquare)) == 1-xsquare
@test log((1-x0)^2) == 2*log(1-x0)
@test real(exp(xI)) == cos(x0)
@test imag(exp(xI)) == sin(x0)
@test exp(xI') == cos(x0)-im*sin(x0)
@test (exp(x0))^(2im) == cos(2x0)+im*sin(2x0)
@test (exp(x0))^Taylor(-5.2im) == cos(5.2x0)-im*sin(5.2x0)
@test abs((tan(x0)).coeffs[8]- 17/315) < tol1
@test abs((tan(x0)).coeffs[14]- 21844/6081075) < tol1
@test evalTaylor(exp(Taylor([0,1],17)),1.0) == e

@test deriv( exp(x(1.0)), 5 ) == exp(1.0)
@test deriv( exp(x(1.0pi)), 3 ) == exp(1.0pi)
@test isapprox( deriv(exp(x(1.0pi)), 10) , exp(1.0pi) )

@test_throws ErrorException 1/x0
@test_throws ErrorException z/z
@test_throws ErrorException x0^1.5
@test_throws ErrorException sqrt(x0)
@test_throws ErrorException log(x0)
@test_throws ErrorException cos(x0)/sin(x0)
@test_throws ErrorException deriv( exp(x(1.0pi)), 30 )

println("    \033[32;1mSUCCESS\033[0m")

