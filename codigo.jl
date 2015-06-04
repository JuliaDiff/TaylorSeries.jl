
using TaylorSeries
Taylor1([1, 2, 3]) # Polynomial of order 2 with coefficients 1, 2, 3
Taylor1([0.0, 1im]) # Also works with complex numbers
affine(a) = a + taylor1_variable(typeof(a),5)  ## a + t of order 5
t = affine(0.0) # Independent variable `t`


t*(3t+2.5)
1/(1-t)
t*(t^2-4)/(t+2)
tI = im*t
t^6  # order is 5
(1-t)^3.2
(1+t)^t
t^3.2


exp(t)
log(1-t)
sqrt(t)
sqrt(1 + t)
imag(exp(tI)')
real(exp(Taylor1([0.0,1im],17))) - cos(Taylor1([0.0,1.0],17)) == 0.0
convert(Taylor1{Rational{Int64}}, exp(t))  # output differes in v0.4


diffTaylor(exp(t))
integTaylor(exp(t))
integTaylor( ans, 1.0)
integTaylor( diffTaylor( exp(-t)), 1.0 ) == exp(-t)
deriv( exp(affine(1.0))) == exp(1.0)
deriv( exp(affine(1.0)), 5) == exp(1.0) # Fifth derivative of `exp(1+t)`


evalTaylor(exp(affine(1.0))) - e # exp(t) around t0=1 (order 5), evaluated there (dt=0)
evalTaylor(exp(t), 1) - e # exp(t) around t0=0 (order 5), evaluated at t=1
evalTaylor( exp( taylor1_variable(17) ), 1) - e # exp(t) around t0=0 (order 17), evaluated at t=1
tBig = Taylor1([zero(BigFloat),one(BigFloat)],50) # With BigFloats
evalTaylor( exp(tBig), one(BigFloat) )
e - ans


show_params_TaylorN()
set_params_TaylorN(8,2)


x = taylorN_variable(1)
y = taylorN_variable(2,6)
typeof(x)
x.order
x.coeffs
get_maxOrder(y)


HomogeneousPolynomial([1,-1])
TaylorN( [HomogeneousPolynomial([1,0]), HomogeneousPolynomial([1,2,3])], 4)


exy = exp(x+y)


get_coeff(exy, [3,5]) == 1/720
rationalize(get_coeff(exy, [3,5]))


f(x,y) = x^3 + 2x^2 * y - 7x + 2
g(x,y) = y - x^4
diffTaylor( f(x,y), 1 )   # partial derivative with respect to 1st variable
diffTaylor( g(x,y), 2 )
diffTaylor( g(x,y), 3 )   # error, since we are dealing with 2 variables


evalTaylor(x+y, [t, 2t])  # x+y, with x=t, y=2t
evalTaylor(exy, [t,2t])
exp(3.0*taylor1_variable(8))


f1 = f(x,y)
g1 = g(x,y)
∇(f1)
gradient( g1 )
jacobian([f1,g1], [2,1])
fg = f1-g1-2*f1*g1
hessian(ans)
evalTaylor(fg, [x+1.0, y+1.0])
hessian(fg, [1.0,1.0])





set_params_TaylorN(4,8)

for i=1:4
    ai = symbol(string("a",i))
    bi = symbol(string("b",i))
    @eval ($ai) = taylorN_variable(Int,$i,4)
    @eval ($bi) = taylorN_variable(Int,4+($i),4)
end
a1
b1


expr_lhs1 = a1^2 + a2^2 + a3^2 + a4^2 ;
expr_lhs2 = b1^2 + b2^2 + b3^2 + b4^2 ;
expr_rhs1 = (a1*b1 - a2*b2 - a3*b3 - a4*b4)^2 ;
expr_rhs2 = (a1*b2 + a2*b1 + a3*b4 - a4*b3)^2 ;
expr_rhs3 = (a1*b3 - a2*b4 + a3*b1 + a4*b2)^2 ;
expr_rhs4 = (a1*b4 + a2*b3 - a3*b2 + a4*b1)^2 ;


lhs = expr_lhs1 * expr_lhs2
rhs = expr_rhs1 + expr_rhs2 + expr_rhs3 + expr_rhs4
lhs == rhs




set_params_TaylorN(40,4)   # maxOrder = 40; numVars = 4
function fateman1(ndeg::Int)
    T = Int128
    unoH = HomogeneousPolynomial(one(T), 0)
    s = TaylorN( [unoH, HomogeneousPolynomial([one(T),one(T),one(T),one(T)],1)], ndeg )
    s = s^ndeg
    # s is converted to order 2*ndeg
    s = TaylorN(s, 2ndeg)
    return s * (s+TaylorN(unoH, 2*ndeg))
end
@time f1 = fateman1(0);
@time f1 = fateman1(20);
get_coeff(f1,[1,6,7,20])
ans > typemax(Int)  # this is the reason for using Int128
function fateman2(ndeg::Int)
    T = Int128
    unoH = HomogeneousPolynomial(one(T), 0)
    s = TaylorN( [unoH, HomogeneousPolynomial([one(T),one(T),one(T),one(T)],1)], ndeg )
    s = s^ndeg
    # s is converted to order 2*ndeg
    s = TaylorN(s, 2ndeg)
    return s^2 + s
end
@time f2 = fateman2(0);
@time f2 = fateman2(20);
get_coeff(f2,[1,6,7,20])
sum(TaylorSeries.sizeTable) # number of distinct monomials
