# Simple tests for TaylorSeries implementation
using TaylorSeries
using Base.Test

x = Taylor([0,1],15)
xI = im*x
z = zero(x)
u = 1.0*one(x)
tol(x) = eps(x)
tol1 = eps(1.0)

@test eltype(convert(Taylor{Complex128},u)) == Complex128
@test eltype(promote(z,Taylor(u))[1]) == Float64
@test eltype(promote(1.0+im, z)[1]) == Complex{Float64}
@test eltype(TaylorSeries.fixshape(z,u)[1]) == Float64
@test length(TaylorSeries.fixshape(Taylor(0,5),z)[1]) == 15
@test TaylorSeries.firstnonzero(x) == 1
@test TaylorSeries.firstnonzero(z) == z.order+1

@test u == 1
@test 0.0 == z
@test x.coeffs[2] == 1
@test z+1 == u
@test x+x == 2x
@test x-x == z

x2 = Taylor([0,0,1],15)
@test x*x == x2
@test (-x)^2 == x2
@test x2/x == x
@test x/(3x) == (1/3)*u
@test x/3im == -xI/3
@test Taylor([0,1,1])/x == x+1
@test (x+im)^2 == x2+2im*x-1
@test imag(x2+2im*x-1) == 2x
@test (Rational(1,2)*x2).coeffs[3] == 1//2
@test x^2/x2 == u
@test ((1+x)^(1/3)).coeffs[3]+1/9 <= tol1
@test 1-x2 == (1+x)-x*(1+x)
@test (1-x2)^2 == (1+x)^2 * (1-x)^2
@test (sqrt(1+x)).coeffs[3] == -1/8
@test ((1-x)^(1/4)).coeffs[15] == -4188908511/549755813888
@test abs(((1+x)^3.2).coeffs[14] + 5.4021062656e-5) < tol1

@test log(exp(x2)) == x2
@test exp(log(1-x2)) == 1-x2
@test log((1-x)^2) == 2*log(1-x)
@test real(exp(xI)) == cos(x)
@test imag(exp(xI)) == sin(x)
@test exp(xI') == cos(x)-im*sin(x)
@test abs((tan(x)).coeffs[8]- 17/315) < tol1
@test abs((tan(x)).coeffs[14]- 21844/6081075) < tol1
@test evalTaylor(exp(Taylor([0,1],17)),1.0) == e

@test_throws t1/x
@test_throws z/z
@test_throws x^1.5
@test_throws sqrt(x)
@test_throws log(x)
@test_throws cos(x)/sin(x)

println("    \033[32;1mSUCCESS\033[0m")

