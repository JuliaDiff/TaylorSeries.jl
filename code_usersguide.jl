
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


evaluate(exp(affine(1.0))) - e # exp(t) around t0=1 (order 5), evaluated there (dt=0)
evaluate(exp(t), 1) - e # exp(t) around t0=0 (order 5), evaluated at t=1
evaluate( exp( taylor1_variable(17) ), 1) - e # exp(t) around t0=0 (order 17), evaluated at t=1
tBig = Taylor1([zero(BigFloat),one(BigFloat)],50) # With BigFloats
evaluate( exp(tBig), one(BigFloat) )
e - ans


x, y = set_variables("x y")


x
typeof(x)
x.order
x.coeffs


set_variables("x y", order=10)


set_variables("α", numvars=3)


show_params_TaylorN()


set_variables("x", numvars=2);
HomogeneousPolynomial([1,-1])
TaylorN( [HomogeneousPolynomial([1,0]), HomogeneousPolynomial([1,2,3])], 4)


x, y = set_variables("x", numvars=2, order=10);
exy = exp(x+y)


get_coeff(exy, [3,5]) == 1/720
rationalize(get_coeff(exy, [3,5]))


f(x,y) = x^3 + 2x^2 * y - 7x + 2
g(x,y) = y - x^4
diffTaylor( f(x,y), 1 )   # partial derivative with respect to 1st variable
diffTaylor( g(x,y), 2 )
diffTaylor( g(x,y), 3 )   # error, since we are dealing with 2 variables



# NEEDS RE-IMPLEMENTATION
# evaluate(x+y, [t, 2t])  # x+y, with x=t, y=2t
# evaluate(exy, [t,2t])
# exp(3.0*taylor1_variable(8))


f1 = f(x,y)
g1 = g(x,y)
∇(f1)
gradient( g1 )
jacobian([f1,g1], [2,1])
fg = f1-g1-2*f1*g1
hessian(ans)
fg1 = f(x+1.0,y+1.0)-g(x+1.0,y+1.0)-2*f(x+1.0,y+1.0)*g(x+1.0,y+1.0)
hessian(fg, [1.0,1.0])
ans == hessian(fg1)





# set_params_TaylorN(4,8)
#
# for i=1:4
#     ai = symbol(string("a",i))
#     bi = symbol(string("b",i))
#     @eval ($ai) = taylorN_variable(Int,$i,4)
#     @eval ($bi) = taylorN_variable(Int,4+($i),4)
# end
# a1
# b1

# Define the variables α₁, ..., α₄, β₁, ..., β₄
make_variable(name, index::Int) = string(name, TaylorSeries.subscriptify(index))
variable_names = [make_variable("α", i) for i in 1:4]
append!(variable_names, [make_variable("β", i) for i in 1:4])
# Create the Taylor objects (order 4, numvars=8)
a1, a2, a3, a4, b1, b2, b3, b4 = set_variables(variable_names, order=4);
a1


lhs1 = a1^2 + a2^2 + a3^2 + a4^2 ;
lhs2 = b1^2 + b2^2 + b3^2 + b4^2 ;
rhs1 = (a1*b1 - a2*b2 - a3*b3 - a4*b4)^2 ;
rhs2 = (a1*b2 + a2*b1 + a3*b4 - a4*b3)^2 ;
rhs3 = (a1*b3 - a2*b4 + a3*b1 + a4*b2)^2 ;
rhs4 = (a1*b4 + a2*b3 - a3*b2 + a4*b1)^2 ;


lhs = lhs1 * lhs2
rhs = rhs1 + rhs2 + rhs3 + rhs4
lhs == rhs




# change number of variables and maxOrder
set_variables("x", numvars=4, order=40)
function fateman1(degree::Int)
    T = Int128
    oneH = HomogeneousPolynomial(one(T), 0)
    # s = 1 + x + y + z + w
    s = TaylorN( [oneH, HomogeneousPolynomial([one(T),one(T),one(T),one(T)],1)], degree )
    s = s^degree
    # s is converted to order 2*ndeg
    s = TaylorN(s, 2*degree)

    s * ( s+TaylorN(oneH, 2*degree) )
end
@time f1 = fateman1(0);
@time f1 = fateman1(20);
get_coeff(f1,[1,6,7,20])
ans > typemax(Int)  # this is the reason for using Int128


function fateman2(degree::Int)
    T = Int128
    oneH = HomogeneousPolynomial(one(T), 0)
    s = TaylorN( [oneH, HomogeneousPolynomial([one(T),one(T),one(T),one(T)],1)], degree )
    s = s^degree
    # s is converted to order 2*ndeg
    s = TaylorN(s, 2*degree)
    return s^2 + s
end
@time f2 = fateman2(0);
@time f2 = fateman2(20);
get_coeff(f2,[1,6,7,20])
sum(TaylorSeries.sizeTable) # number of distinct monomials
