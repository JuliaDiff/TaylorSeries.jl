var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#TaylorSeries.jl-1",
    "page": "Home",
    "title": "TaylorSeries.jl",
    "category": "section",
    "text": "A Julia package for Taylor expansions in one or more independent variables."
},

{
    "location": "index.html#Authors-1",
    "page": "Home",
    "title": "Authors",
    "category": "section",
    "text": "Luis Benet, Instituto de Ciencias   Físicas, Universidad Nacional Autónoma de México (UNAM).\nDavid P. Sanders, Facultad   de Ciencias, Universidad Nacional Autónoma de México (UNAM)."
},

{
    "location": "index.html#License-1",
    "page": "Home",
    "title": "License",
    "category": "section",
    "text": "TaylorSeries is licensed under the MIT \"Expat\" license; see LICENSE for the full license text."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "TaylorSeries.jl is a registered package, and is simply installed by runningjulia> Pkg.add(\"TaylorSeries\")"
},

{
    "location": "index.html#Related-packages-1",
    "page": "Home",
    "title": "Related packages",
    "category": "section",
    "text": "Polynomials.jl: Polynomial manipulations\nPowerSeries.jl: Truncated power series for Julia\nMultiPoly.jl Sparse multivariate polynomials in Julia"
},

{
    "location": "index.html#Acknowledgments-1",
    "page": "Home",
    "title": "Acknowledgments",
    "category": "section",
    "text": "This project began (using Python) during a Masters\' course in the postgraduate programs in Physics and in Mathematics at UNAM, during the second half of 2013. We thank the participants of the course for putting up with the half-baked material and contributing energy and ideas.We acknowledge financial support from DGAPA-UNAM PAPIME grants PE-105911 and PE-107114, and PAPIIT grants IG-101113 and IG-100616. LB acknowledges support through a Cátedra Marcos Moshinsky (2013)."
},

{
    "location": "background.html#",
    "page": "Background",
    "title": "Background",
    "category": "page",
    "text": ""
},

{
    "location": "background.html#Background-1",
    "page": "Background",
    "title": "Background",
    "category": "section",
    "text": ""
},

{
    "location": "background.html#Introduction-1",
    "page": "Background",
    "title": "Introduction",
    "category": "section",
    "text": "TaylorSeries.jl is an implementation of high-order automatic differentiation, as presented in the book by W. Tucker [1]. The general idea is the following.The Taylor series expansion of an analytical function f(t) with one independent variable t around t_0 can be written asbeginequation\nf(t) = f_0 + f_1 (t-t_0) + f_2 (t-t_0)^2 + cdots + f_k (t-t_0)^k + cdots\nendequationwhere f_0=f(t_0), and the Taylor coefficients f_k = f_k(t_0) are the k-th normalized derivatives at t_0:beginequation\nf_k = frac1k fracrm d^k f rm d t^k(t_0)\nendequationThus, computing the high-order derivatives of f(t) is equivalent to computing its Taylor expansion.In the case of many independent variables the same statements hold, though things become more subtle. Following Alex Haro\'s approach [2], the Taylor expansion is an infinite sum of homogeneous polynomials in the d independent variables x_1 x_2 dots x_d, which takes the formbeginequation\nf_k (mathbfx_0) = sum_m_1+cdots+m_d = k f_m_1dotsm_d \n(x_1-x_0_1)^m_1 cdots (x_d-x_0_d)^m_d =\nsum_mathbfm=k f_mathbfm (mathbfx-mathbfx_0)^mathbfm\nendequationHere, mathbfmin mathbbN^d is a multi-index of the k-th order homogeneous polynomial and mathbfx=(x_1x_2ldotsx_d) are the d independent variables.In both cases, a Taylor series expansion can be represented by a vector containing its coefficients. The difference between the cases of one or more independent variables is that the coefficients are real or complex numbers in the former case, but homogeneous polynomials in the latter case. This motivates the construction of the Taylor1 and TaylorN types."
},

{
    "location": "background.html#Arithmetic-operations-1",
    "page": "Background",
    "title": "Arithmetic operations",
    "category": "section",
    "text": "Arithmetic operations involving Taylor series can be expressed as operations on the coefficients:begineqnarray\nlabeleqarith1\n(f(x) pm g(x))_k  =  f_k pm g_k  \nlabeleqarith2\n(f(x) cdot g(x))_k  =  sum_i=0^k f_i  g_k-i  \nlabeleqarith3\nBig( fracf(x)g(x) Big)_k  =  frac1g_0 Big f_k -\nsum_i=0^k-1 big(fracf(x)g(x)big)_i  g_k-i Big \nendeqnarrayEquations (\\ref{eq:arith1}-\\ref{eq:arith3}) corresponds to a convolution."
},

{
    "location": "background.html#Elementary-functions-of-polynomials-1",
    "page": "Background",
    "title": "Elementary functions of polynomials",
    "category": "section",
    "text": "Consider a function y(t) that satisfies the ordinary differential equation doty = f(y), y(t_0)=y_0, where t is the independent variable. Writing y(t) and f(t) as Taylor polynomials of t, substituting these in the differential equation and equating equal powers of the independent variable leads to the recursion relationbeginequation\nlabeleqrec\ny_n+1 = fracf_nn+1\nendequationEquation (\\ref{eq:rec}) and the corresponding initial condition y(t_0)=y_0 define a recurrence relation for the Taylor coefficients of y(t) around t_0.The following are  examples of such recurrence relations for some elementary functions:begineqnarray\np(t)=(f(t))^alpha  qquad \n  p_k  = frac1k  f_0sum_j=0^k-1big( alpha(k-j)-jbig)\n   f_k-j  p_j \ne(t) = exp(f(t))  qquad \n  e_k  = frac1ksum_j=0^k-1 (k-j)  f_k-j  e_j \nl(t) = log(f(t))  qquad \n  l_k  = frac1f_0big( f_k - frac1ksum_j=1^k-1 j\n     f_k-j  l_j big) \ns(t) = sin(f(t))  qquad \n  s_k  = frac1ksum_j=0^k-1 (k-j)  f_k-j  c_j \nc(t) = cos(f(t))  qquad \n  c_k  = -frac1ksum_j=0^k-1 (k-j)  f_k-j  s_j\nendeqnarrayThe recursion relations for s(t) = sinbig(f(t)big) and c(t) = cosbig(f(t)big) depend on each other; this reflects the fact that they are solutions of a second-order differential equation.All these relations hold for Taylor expansions in one and more independent variables; in the latter case, the Taylor coefficients f_k are homogeneous polynomials of degree k; see [2]."
},

{
    "location": "background.html#refs-1",
    "page": "Background",
    "title": "References",
    "category": "section",
    "text": "[1] W. Tucker, Validated Numerics: A Short Introduction to Rigorous Computations, Princeton University Press (2011).[2] A. Haro, Automatic differentiation methods in computational dynamical systems: Invariant manifolds and normal forms of vector fields at fixed points, preprint."
},

{
    "location": "userguide.html#",
    "page": "User guide",
    "title": "User guide",
    "category": "page",
    "text": ""
},

{
    "location": "userguide.html#User-guide-1",
    "page": "User guide",
    "title": "User guide",
    "category": "section",
    "text": "CurrentModule = TaylorSeriesTaylorSeries.jl is a basic polynomial algebraic manipulator in one or more variables; these two cases are treated separately.  Three new types are defined, Taylor1, HomogeneousPolynomial and TaylorN, which correspond to expansions in one independent variable, homogeneous polynomials of various variables, and the polynomial series in many independent variables, respectively. These types are subtypes of AbstractSeries, which in turn is a subtype of Number, and are defined parametrically.The package is loaded as usual:using TaylorSeries"
},

{
    "location": "userguide.html#One-independent-variable-1",
    "page": "User guide",
    "title": "One independent variable",
    "category": "section",
    "text": "Taylor expansions in one variable are represented by the Taylor1 type, which consists of a vector of coefficients (fieldname coeffs) and the maximum order considered for the expansion (fieldname order). The coefficients are arranged in ascending order with respect to the degree of the monomial, so that coeffs[1] is the constant term, coeffs[2] gives the first order term (t^1), etc. Yet, it is possible to have the natural ordering with respect to the degree; see below. This is a dense representation of the polynomial. The order of the polynomial can be omitted in the constructor, which is then fixed by the length of the vector of coefficients. If the length of the vector does not correspond with the order, order is used, which effectively truncates polynomial to degree order.Taylor1([1, 2, 3],4) # Polynomial of order 4 with coefficients 1, 2, 3\nTaylor1([0.0, 1im]) # Also works with complex numbers\nTaylor1(ones(8), 2) # Polynomial truncated to order 2\nshift_taylor(a) = a + Taylor1(typeof(a),5)  ## a + taylor-polynomial of order 5\nt = shift_taylor(0.0) # Independent variable `t`Note that the information about the maximum order considered is displayed using a big-𝒪 notation. In some cases, it is desirable to not display the big-𝒪 notation. The function displayBigO allows to control whether it is displayed or not.displayBigO(false) # turn-off displaying big O notation\nt\ndisplayBigO(true) # turn it on\ntSimilarly, it is possible to control if the format of the displayed series through the function use_show_default; use_show_default(true) uses the Base.show_default, while use_show_default(false) uses the custom display form (default).use_show_default(true) # use Base.show method\nt\nuse_show_default(false) # use custom `show`\ntThe definition of shift_taylor(a) uses the method Taylor1([::Type{Float64}], [order::Int64=1]), which is a shortcut to define the independent variable of a Taylor expansion, of given type and order (defaults are Float64 and order=1). This is one of the easiest ways to work with the package.The usual arithmetic operators (+, -, *, /, ^, ==) have been extended to work with the Taylor1 type, including promotions that involve Numbers. The operations return a valid Taylor expansion of maximum order. This is apparent in the last example below, where the answer is beyond the order of the expansion.t*(3t+2.5)\n1/(1-t)\nt*(t^2-4)/(t+2)\ntI = im*t\n(1-t)^3.2\n(1+t)^t\nt^6  # t is of order 5If no valid Taylor expansion can be computed, an error is thrown, for instance when a derivative is not defined (or simply diverges):1/t\nt^3.2\nabs(t)Several elementary functions have been implemented; their coefficients are computed recursively. At the moment of this writing, these functions are exp, log, sqrt, the trigonometric functions sin, cos and tan, their inverses, as well as the hyperbolic functions sinh, cosh and tanh and their inverses; more will be added in the future. Note that this way of obtaining the Taylor coefficients is not the laziest way, in particular for many independent variables. Yet, it is quite efficient, especially for the integration of ordinary differential equations, which is among the applications we have in mind (see also TaylorIntegration.jl).exp(t)\nlog(1-t)\nsqrt(1 + t)\nimag(exp(tI)\')\nreal(exp(Taylor1([0.0,1im],17))) - cos(Taylor1([0.0,1.0],17)) == 0.0\nconvert(Taylor1{Rational{Int64}}, exp(t))Again, errors are thrown whenever it is necessary.sqrt(t)\nlog(t)To obtain a specific coefficient, getcoeff can be used. Another alternative is to request the specific degree using the vector notation, where the index corresponds to the degree of the term.expon = exp(t)\ngetcoeff(expon, 0) == expon[0]\nrationalize(expon[3])Differentiating and integrating is straightforward for polynomial expansions in one variable, using derivative and integrate. These functions return the corresponding Taylor1 expansions. The last coefficient of a derivative is set to zero to keep the same order as the original polynomial; for the integral, an integration constant may be set by the user (the default is zero). The order of the resulting polynomial is not changed. The value of the n-th (n ge 0) derivative is obtained using derivative(n,a), where a is a Taylor series.derivative(exp(t))\nintegrate(exp(t))\nintegrate( exp(t), 1.0)\nintegrate( derivative( exp(-t)), 1.0 ) == exp(-t)\nderivative(1, exp(shift_taylor(1.0))) == exp(1.0)\nderivative(5, exp(shift_taylor(1.0))) == exp(1.0) # 5-th derivative of `exp(1+t)`To evaluate a Taylor series at a given point, Horner\'s rule is used via the function evaluate(a, dt). Here, dt is the increment from the point t_0 around which the Taylor expansion of a is calculated, i.e., the series is evaluated at t = t_0 + dt. Omitting dt corresponds to dt = 0; see evaluate.evaluate(exp(shift_taylor(1.0))) - ℯ # exp(t) around t0=1 (order 5), evaluated there (dt=0)\nevaluate(exp(t), 1) - ℯ # exp(t) around t0=0 (order 5), evaluated at t=1\nevaluate(exp( Taylor1(17) ), 1) - ℯ # exp(t) around t0=0, order 17\ntBig = Taylor1(BigFloat, 50) # Independent variable with BigFloats, order 50\neBig = evaluate( exp(tBig), one(BigFloat) )\nℯ - eBigAnother way to obtain the value of a Taylor1 polynomial p at a given value x, is to call p as if it were a function, i.e., p(x):t = Taylor1(15)\np = sin(t)\nevaluate(p, pi/2) # value of p at pi/2 using `evaluate`\np(pi/2) # value of p at pi/2 by evaluating p as a function\np(pi/2) == evaluate(p, pi/2)\np(0.0)\np() == p(0.0) # p() is a shortcut to obtain the 0-th order coefficient of `p`Note that the syntax p(x) is equivalent to evaluate(p, x), whereas p() is equivalent to evaluate(p). For more details about function-like behavior for a given type in Julia, see the Function-like objects section of the Julia manual.Useful shortcuts are taylor_expand are update!. The former returns the expansion of a function around a given value t0. In turn, update! provides an in-place update of a given Taylor polynomial, that is, it shifts it further by the provided amount.p = taylor_expand( x -> sin(x), pi/2, order=16) # 16-th order expansion of sin(t) around pi/2\nupdate!(p, 0.025) # updates the expansion given by p, by shifting it further by 0.025\np"
},

{
    "location": "userguide.html#Many-variables-1",
    "page": "User guide",
    "title": "Many variables",
    "category": "section",
    "text": "A polynomial in N1 variables can be represented in (at least) two ways: As a vector whose coefficients are homogeneous polynomials of fixed degree, or as a vector whose coefficients are polynomials in N-1 variables. The current implementation of TaylorSeries.jl corresponds to the first option, though some infrastructure has been built that permits to develop the second one. An elegant (lazy) implementation of the second representation was discussed  here.The structure TaylorN is constructed as a vector of parameterized homogeneous polynomials defined by the type HomogeneousPolynomial, which in turn is a vector of coefficients of given order (degree). This implementation imposes the user to specify the (maximum) order considered and the number of independent variables at the beginning, which can be conveniently done using set_variables. A vector of the resulting Taylor variables is returned:x, y = set_variables(\"x y\")\ntypeof(x)\nx.order\nx.coeffsAs shown, the resulting objects are of TaylorN{Float64} type. There is an optional order keyword argument in set_variables, used to specify the maximum order of the TaylorN polynomials. Note that one can specify the variables using a vector of symbols.set_variables([:x, :y], order=10)Similarly, numbered variables are also available by specifying a single variable name and the optional keyword argument numvars:set_variables(\"α\", numvars=3)Alternatively to set_variables, get_variables can be used if one doesn\'t want to change internal dictionaries. get_variables returns a vector of TaylorN independent variables of a desired order (lesser than get_order so the internals doesn\'t have to change) with the length and variable names defined by set_variables initially.get_variables(2) # vector of independent variables of order 2The function show_params_TaylorN displays the current values of the parameters, in an info block.show_params_TaylorN()Internally, changing the parameters (maximum order and number of variables) redefines the hash-tables that translate the index of the coefficients of a HomogeneousPolynomial of given order into the corresponding multi-variable monomials, or the other way around. Fixing these values from the start is imperative; the initial (default) values are order = 6 and num_vars=2.The easiest way to construct a TaylorN object is by defining the independent variables. This can be done using set_variables as above, or through the method TaylorN{T<:Number}(::Type{T}, nv::Int) for the nv independent TaylorN{T} variable; the order can be also specified using the optional keyword argument order.x, y = set_variables(\"x y\", numvars=2, order=6);\nx\nTaylorN(1, order=4) # variable 1 of order 4\nTaylorN(Int, 2)    # variable 2, type Int, order=get_order()=6Other ways of constructing TaylorN polynomials involve using HomogeneousPolynomial objects directly, which is uncomfortable.set_variables(:x, numvars=2); # symbols can be used\nHomogeneousPolynomial([1,-1])\nTaylorN([HomogeneousPolynomial([1,0]), HomogeneousPolynomial([1,2,3])],4)The Taylor expansions are implemented around 0 for all variables; if the expansion is needed around a different value, the trick is a simple translation of the corresponding independent variable, i.e. x to x+a.As before, the usual arithmetic operators (+, -, *, /, ^, ==) have been extended to work with TaylorN objects, including the appropriate promotions to deal with numbers. (Some of the arithmetic operations have been extended for HomogeneousPolynomial, whenever the result is a HomogeneousPolynomial; division, for instance, is not extended.)Also, the elementary functions have been implemented, again by computing their coefficients recursively:x, y = set_variables(\"x y\", order=10);\nexy = exp(x+y)The function getcoeff gives the normalized coefficient of the polynomial that corresponds to the monomial specified by the tuple or vector v containing the powers. For instance, for the polynomial exy above, the coefficient of the monomial x^3 y^5 is obtained using getcoeff(exy, (3,5)) or getcoeff(exy, [3,5]).getcoeff(exy, (3,5))\nrationalize(ans)Similar to Taylor1, vector notation can be used to request specific coefficients of HomogeneousPolynomial or TaylorN objects. For TaylorN objects, the index refers to the degree of the HomogeneousPolynomial. In the case of HomogeneousPolynomial the index refers to the position of the hash table. The function show_monomials can be used to obtain the coefficient a specific monomial, given the degree of the HomogeneousPolynomial.exy[8] # get the 8th order term\nshow_monomials(8)\nexy[8][6] # get the 6th coeff of the 8th order termPartial differentiation is also implemented for TaylorN objects, through the function derivative, specifying the number of the variable, or its symbol, as the second argument.p = x^3 + 2x^2 * y - 7x + 2\nq = y - x^4\nderivative( p, 1 )   # partial derivative with respect to 1st variable\nderivative( q, :y )  # partial derivative with respect to :yIf we ask for the partial derivative with respect to a non-defined variable, an error is thrown.derivative( q, 3 )   # error, since we are dealing with 2 variablesTo obtain more specific partial derivatives we have two specialized methods that involve a tuple, which represents the number of derivatives with respect to each variable (so the tuple\'s length has to be the same as the actual number of variables). These methods either return the TaylorN object in question, or the coefficient corresponding to the specified tuple, normalized by the factorials defined by the tuple. The latter is in essence the 0-th order coefficient of the former.derivative(p, (2,1)) # two derivatives on :x and one on :y\nderivative((2,1), p) # 0-th order coefficient of the previous expression\nderivative(p, (1,1)) # one derivative on :x and one on :y\nderivative((1,1), p) # 0-th order coefficient of the previous expressionIntegration with respect to the r-th variable for HomogeneousPolynomials and TaylorN objects is obtained using integrate. Note that integrate for TaylorN objects allows to specify a constant of integration, which must be independent from the integrated variable. Again, the integration variable may be specified by its symbol.integrate( derivative( p, 1 ), 1) # integrate with respect to the first variable\nintegrate( derivative( p, 1 ), :x, 2) # integration with respect to :x, constant of integration is 2\nintegrate( derivative( q, 2 ), :y, -x^4) == q\nintegrate( derivative( q, 2 ), 2, y)evaluate can also be used for TaylorN objects, using it on vectors of numbers (Real or Complex); the length of the vector must coincide with the number of independent variables. evaluate also allows to specify only one variable and a value.evaluate(exy, [.1,.02]) == exp(0.12)\nevaluate(exy, :x, 0.0) == exp(y)  # evaluate `exy` for :x -> 0Analogously to Taylor1, another way to obtain the value of a TaylorN polynomial p at a given point x, is to call it as if it were a function: the syntax p(x) for p::TaylorN is equivalent to evaluate(p,x), and p() is equivalent to evaluate(p).exy([.1,.02]) == exp(0.12)\nexy(:x, 0.0)The functions taylor_expand and update! work as well for TaylorN.xysq = x^2 + y^2\nupdate!(xysq, [1.0, -2.0]) # expand around (1,-2)\nxysq\nupdate!(xysq, [-1.0, 2.0]) # shift-back\nxysq == x^2 + y^2Functions to compute the gradient, Jacobian and Hessian have also been implemented. Using the polynomials p = x^3 + 2x^2 y - 7 x + 2 and q = y-x^4 defined above, we may use gradient (or ∇); the results are of type Array{TaylorN{T},1}. To compute the Jacobian and Hessian of a vector field evaluated at a point, we use respectively jacobian and hessian:∇(p)\ngradient( q )\nr = p-q-2*p*q\nhessian(ans)\njacobian([p,q], [2,1])\nhessian(r, [1.0,1.0])Other specific applications are described in the Examples."
},

{
    "location": "userguide.html#Mixtures-1",
    "page": "User guide",
    "title": "Mixtures",
    "category": "section",
    "text": "As mentioned above, Taylor1{T}, HomogeneousPolynomial{T} and TaylorN{T} are parameterized structures such that T<:AbstractSeries, the latter is a subtype of Number. Then, we may actually define Taylor expansions in N+1 variables, where one of the variables (the Taylor1 variable) is somewhat special.x, y = set_variables(\"x y\", order=3)\nt1N = Taylor1([zero(x), one(x)], 5)The last line defines a Taylor1{TaylorN{Float64}} variable, which is of order 5 in t and order 3 in x and y. Then, we can evaluate functions involving such polynomials:cos(2.1+x+t1N)This kind of expansions are of interest when studying the dependence of parameters, for instance in the context of bifurcation theory or when considering the dependence of the solution of a differential equation on the initial conditions, around a given solution. In this case, x and y represent small variations around a given value of the parameters, or around some specific initial condition. Such constructions are exploited in the package TaylorIntegration.jl.CurrentModule = nothing"
},

{
    "location": "examples.html#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "examples.html#Examples-1",
    "page": "Examples",
    "title": "Examples",
    "category": "section",
    "text": "CurrentModule = TaylorSeries"
},

{
    "location": "examples.html#Four-square-identity-1",
    "page": "Examples",
    "title": "Four-square identity",
    "category": "section",
    "text": "The first example shows that the four-square identity holds:begineqnarray\n(a_1+a_2+a_3+a_4)cdot(b_1+b_2+b_3+b_4)  = \n     (a_1 b_1 - a_2 b_2 - a_3 b_3 -a_4 b_4)^2 + qquad nonumber \nlabeleqEuler\n    (a_1 b_2 - a_2 b_1 - a_3 b_4 -a_4 b_3)^2 + \n    (a_1 b_3 - a_2 b_4 - a_3 b_1 -a_4 b_2)^2 + nonumber \n    (a_1 b_4 - a_2 b_3 - a_3 b_2 -a_4 b_1)^2 nonumber\nendeqnarraywhich was originally proved by Euler. The code can also be found in this test of the package.First, we reset the maximum degree of the polynomial to 4, since the RHS of the equation has a priori terms of fourth order, and define the 8 independent variables.using TaylorSeries\n# Define the variables α₁, ..., α₄, β₁, ..., β₄\nmake_variable(name, index::Int) = string(name, TaylorSeries.subscriptify(index))\nvariable_names = [make_variable(\"α\", i) for i in 1:4]\nappend!(variable_names, [make_variable(\"β\", i) for i in 1:4])\n# Create the TaylorN variables (order=4, numvars=8)\na1, a2, a3, a4, b1, b2, b3, b4 = set_variables(variable_names, order=4)\na1 # variable a1Now we compute each term appearing in Eq. (\\ref{eq:Euler})# left-hand side\nlhs1 = a1^2 + a2^2 + a3^2 + a4^2 ;\nlhs2 = b1^2 + b2^2 + b3^2 + b4^2 ;\nlhs = lhs1 * lhs2;\n# right-hand side\nrhs1 = (a1*b1 - a2*b2 - a3*b3 - a4*b4)^2 ;\nrhs2 = (a1*b2 + a2*b1 + a3*b4 - a4*b3)^2 ;\nrhs3 = (a1*b3 - a2*b4 + a3*b1 + a4*b2)^2 ;\nrhs4 = (a1*b4 + a2*b3 - a3*b2 + a4*b1)^2 ;\nrhs = rhs1 + rhs2 + rhs3 + rhs4;We now compare the two sides of the identity,lhs == rhsThe identity is satisfied. square"
},

{
    "location": "examples.html#Fateman-test-1",
    "page": "Examples",
    "title": "Fateman test",
    "category": "section",
    "text": "Richard J. Fateman, from Berkley, proposed as a stringent test of polynomial multiplication the evaluation of s*(s+1), where s = (1+x+y+z+w)^20. This is implemented in the function fateman1 below. We shall also consider the form s^2+s in fateman2, which involves fewer operations (and makes a fairer comparison to what Mathematica does).using TaylorSeries\nconst order = 20\nconst x, y, z, w = set_variables(Int128, \"x\", numvars=4, order=2order)\nfunction fateman1(degree::Int)\n    T = Int128\n    s = one(T) + x + y + z + w\n    s = s^degree\n    s * ( s + one(T) )\nend(In the following lines, which are run when the documentation is built, by some reason the timing appears before the command executed.)@time fateman1(0);\n# hide\n@time f1 = fateman1(20);Another implementation of the same, but exploiting optimizations related to ^2 yields:function fateman2(degree::Int)\n    T = Int128\n    s = one(T) + x + y + z + w\n    s = s^degree\n    s^2 + s\nend\nfateman2(0);\n@time f2 = fateman2(20); # the timing appears aboveWe note that the above functions use expansions in Int128. This is actually required, since some coefficients are larger than typemax(Int):getcoeff(f2, (1,6,7,20)) # coefficient of x y^6 z^7 w^{20}\nans > typemax(Int)\nlength(f2)\nsum(TaylorSeries.size_table)\nset_variables(\"x\", numvars=2, order=6) # hideThese examples show that fateman2 is nearly twice as fast as fateman1, and that the series has 135751 monomials in 4 variables."
},

{
    "location": "examples.html#Bechmarks-1",
    "page": "Examples",
    "title": "Bechmarks",
    "category": "section",
    "text": "The functions described above have been compared against Mathematica v11.1. The relevant files used for benchmarking can be found here. Running on a MacPro with Intel-Xeon processors 2.7GHz, we obtain that Mathematica requires on average (5 runs) 3.075957 seconds for the computation, while for fateman1 and fateman2 above we obtain 2.15408 and 1.08337, respectively.Then, with the current version of TaylorSeries.jl and using Julia v0.7.0, our implementation of fateman1 is about 30%-40% faster. (The original test by Fateman corresponds to fateman1 above, which avoids some optimizations related to squaring; the implementation in Mathematica is done such that this optimization does not occur.)"
},

{
    "location": "api.html#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api.html#Library-1",
    "page": "API",
    "title": "Library",
    "category": "section",
    "text": "CurrentModule = TaylorSeries"
},

{
    "location": "api.html#TaylorSeries.TaylorSeries",
    "page": "API",
    "title": "TaylorSeries.TaylorSeries",
    "category": "module",
    "text": "TaylorSeries\n\nA Julia package for Taylor expansions in one or more independent variables.\n\nThe basic constructors are Taylor1 and TaylorN; see also HomogeneousPolynomial.\n\n\n\n\n\n"
},

{
    "location": "api.html#Module-1",
    "page": "API",
    "title": "Module",
    "category": "section",
    "text": "TaylorSeries"
},

{
    "location": "api.html#TaylorSeries.Taylor1",
    "page": "API",
    "title": "TaylorSeries.Taylor1",
    "category": "type",
    "text": "Taylor1{T<:Number} <: AbstractSeries{T}\n\nDataType for polynomial expansions in one independent variable.\n\nFields:\n\ncoeffs :: Array{T,1} Expansion coefficients; the i-th   component is the coefficient of degree i-1 of the expansion.\norder  :: Int64 Maximum order (degree) of the polynomial.\n\nNote that Taylor1 variables are callable. For more information, see evaluate.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.HomogeneousPolynomial",
    "page": "API",
    "title": "TaylorSeries.HomogeneousPolynomial",
    "category": "type",
    "text": "HomogeneousPolynomial{T<:Number} <: AbstractSeries{T}\n\nDataType for homogenous polynomials in many (>1) independent variables.\n\nFields:\n\ncoeffs  :: Array{T,1} Expansion coefficients of the homogeneous\n\npolynomial; the i-th component is related to a monomial, where the degrees of the independent variables are specified by coeff_table[order+1][i].\n\norder   :: Int order (degree) of the homogenous polynomial.\n\nNote that HomogeneousPolynomial variables are callable. For more information, see evaluate.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.TaylorN",
    "page": "API",
    "title": "TaylorSeries.TaylorN",
    "category": "type",
    "text": "TaylorN{T<:Number} <: AbstractSeries{T}\n\nDataType for polynomial expansions in many (>1) independent variables.\n\nFields:\n\ncoeffs  :: Array{HomogeneousPolynomial{T},1} Vector containing the\n\nHomogeneousPolynomial entries. The i-th component corresponds to the homogeneous polynomial of degree i-1.\n\norder   :: Int  maximum order of the polynomial expansion.\n\nNote that TaylorN variables are callable. For more information, see evaluate.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.AbstractSeries",
    "page": "API",
    "title": "TaylorSeries.AbstractSeries",
    "category": "type",
    "text": "AbstractSeries{T<:Number} <: Number\n\nParameterized abstract type for Taylor1, HomogeneousPolynomial and TaylorN.\n\n\n\n\n\n"
},

{
    "location": "api.html#Types-1",
    "page": "API",
    "title": "Types",
    "category": "section",
    "text": "Taylor1\nHomogeneousPolynomial\nTaylorN\nAbstractSeries"
},

{
    "location": "api.html#TaylorSeries.Taylor1-Union{Tuple{Type{T}}, Tuple{T}, Tuple{Type{T},Int64}} where T<:Number",
    "page": "API",
    "title": "TaylorSeries.Taylor1",
    "category": "method",
    "text": "Taylor1([T::Type=Float64], [order::Int=1])\n\nShortcut to define the independent variable of a Taylor1{T} polynomial of given order. The default type for T is Float64.\n\njulia> Taylor1(16)\n 1.0 t + 𝒪(t¹⁷)\n\njulia> Taylor1(Rational{Int}, 4)\n 1//1 t + 𝒪(t⁵)\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.HomogeneousPolynomial-Union{Tuple{T}, Tuple{Type{T},Int64}} where T<:Number",
    "page": "API",
    "title": "TaylorSeries.HomogeneousPolynomial",
    "category": "method",
    "text": "HomogeneousPolynomial([T::Type=Float64], nv::Int])\n\nShortcut to define the nv-th independent HomogeneousPolynomial{T}. The default type for T is Float64.\n\njulia> HomogeneousPolynomial(1)\n 1.0 x₁\n\njulia> HomogeneousPolynomial(Rational{Int}, 2)\n 1//1 x₂\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.TaylorN-Union{Tuple{T}, Tuple{Type{T},Int64}} where T<:Number",
    "page": "API",
    "title": "TaylorSeries.TaylorN",
    "category": "method",
    "text": "TaylorN([T::Type=Float64], nv::Int; [order::Int=get_order()])\n\nShortcut to define the nv-th independent TaylorN{T} variable as a polynomial. The order is defined through the keyword parameter order, whose default corresponds to get_order(). The default of type for T is Float64.\n\njulia> TaylorN(1)\n 1.0 x₁ + 𝒪(‖x‖⁷)\n\njulia> TaylorN(Rational{Int},2)\n 1//1 x₂ + 𝒪(‖x‖⁷)\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.set_variables",
    "page": "API",
    "title": "TaylorSeries.set_variables",
    "category": "function",
    "text": "set_variables([T::Type], names::String; [order=get_order(), numvars=-1])\n\nReturn a TaylorN{T} vector with each entry representing an independent variable. names defines the output for each variable (separated by a space). The default type T is Float64, and the default for order is the one defined globally. Changing the order or numvars resets the hash_tables.\n\nIf numvars is not specified, it is inferred from names. If only one variable name is defined and numvars>1, it uses this name with subscripts for the different variables.\n\njulia> set_variables(Int, \"x y z\", order=4)\n3-element Array{TaylorSeries.TaylorN{Int64},1}:\n  1 x + 𝒪(‖x‖⁵)\n  1 y + 𝒪(‖x‖⁵)\n  1 z + 𝒪(‖x‖⁵)\n\njulia> set_variables(\"α\", numvars=2)\n2-element Array{TaylorSeries.TaylorN{Float64},1}:\n  1.0 α₁ + 𝒪(‖x‖⁵)\n  1.0 α₂ + 𝒪(‖x‖⁵)\n\njulia> set_variables(\"x\", order=6, numvars=2)\n2-element Array{TaylorSeries.TaylorN{Float64},1}:\n  1.0 x₁ + 𝒪(‖x‖⁷)\n  1.0 x₂ + 𝒪(‖x‖⁷)\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.get_variables",
    "page": "API",
    "title": "TaylorSeries.get_variables",
    "category": "function",
    "text": "get_variables(T::Type, [order::Int=get_order()])\n\nReturn a TaylorN{T} vector with each entry representing an independent variable. It takes the default _params_TaylorN_ values if set_variables hasn\'t been changed with the exception that order can be explicitely established by the user without changing internal values for num_vars or variable_names. Ommiting T defaults to Float64.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.show_params_TaylorN",
    "page": "API",
    "title": "TaylorSeries.show_params_TaylorN",
    "category": "function",
    "text": "show_params_TaylorN()\n\nDisplay the current parameters for TaylorN and HomogeneousPolynomial types.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.show_monomials",
    "page": "API",
    "title": "TaylorSeries.show_monomials",
    "category": "function",
    "text": "show_monomials(ord::Int) --> nothing\n\nList the indices and corresponding of a HomogeneousPolynomial of degree ord.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.getcoeff",
    "page": "API",
    "title": "TaylorSeries.getcoeff",
    "category": "function",
    "text": "getcoeff(a, n)\n\nReturn the coefficient of order n::Int of a a::Taylor1 polynomial.\n\n\n\n\n\ngetcoeff(a, v)\n\nReturn the coefficient of a::HomogeneousPolynomial, specified by v, which is a tuple (or vector) with the indices of the specific monomial.\n\n\n\n\n\ngetcoeff(a, v)\n\nReturn the coefficient of a::TaylorN, specified by v, which is a tuple (or vector) with the indices of the specific monomial.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.evaluate",
    "page": "API",
    "title": "TaylorSeries.evaluate",
    "category": "function",
    "text": "evaluate(a, [dx])\n\nEvaluate a Taylor1 polynomial using Horner\'s rule (hand coded). If dx is ommitted, its value is considered as zero. Note that the syntax a(dx) is equivalent to evaluate(a,dx), and a() is equivalent to evaluate(a).\n\n\n\n\n\nevaluate(x, δt)\n\nEvaluates each element of x::Union{ Vector{Taylor1{T}}, Matrix{Taylor1{T}} }, representing the dependent variables of an ODE, at time δt. Note that the syntax x(δt) is equivalent to evaluate(x, δt), and x() is equivalent to evaluate(x).\n\n\n\n\n\nevaluate(a, x)\n\nSubstitute x::Taylor1 as independent variable in a a::Taylor1 polynomial. Note that the syntax a(x) is equivalent to evaluate(a, x).\n\n\n\n\n\nevaluate(a, [vals])\n\nEvaluate a HomogeneousPolynomial polynomial at vals. If vals is ommitted, it\'s evaluated at zero. Note that the syntax a(vals) is equivalent to evaluate(a, vals); and a() is equivalent to evaluate(a).\n\n\n\n\n\nevaluate(a, [vals])\n\nEvaluate the TaylorN polynomial a at vals. If vals is ommitted, it\'s evaluated at zero. Note that the syntax a(vals) is equivalent to evaluate(a, vals); and a() is equivalent to evaluate(a).\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.evaluate!",
    "page": "API",
    "title": "TaylorSeries.evaluate!",
    "category": "function",
    "text": "evaluate!(x, δt, x0)\n\nEvaluates each element of x::Array{Taylor1{T},1}, representing the Taylor expansion for the dependent variables of an ODE at time δt. It updates the vector x0 with the computed values.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.taylor_expand",
    "page": "API",
    "title": "TaylorSeries.taylor_expand",
    "category": "function",
    "text": "taylor_expand(f, x0; order)\n\nComputes the Taylor expansion of the function f around the point x0.\n\nIf x0 is a scalar, a Taylor1 expansion will be returned. If x0 is a vector, a TaylorN expansion will be computed. If the dimension of x0 (length(x0)) is different from the variables set for TaylorN (get_numvars()), an AssertionError will be thrown.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.update!",
    "page": "API",
    "title": "TaylorSeries.update!",
    "category": "function",
    "text": "update!(a, x0)\n\nTakes a <: Union{Taylo1,TaylorN} and expands it around the coordinate x0.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.derivative",
    "page": "API",
    "title": "TaylorSeries.derivative",
    "category": "function",
    "text": "derivative(a)\n\nReturn the Taylor1 polynomial of the differential of a::Taylor1. The last coefficient is set to zero.\n\n\n\n\n\nderivative(a, n)\n\nCompute recursively the Taylor1 polynomial of the n-th derivative of a::Taylor1.\n\n\n\n\n\nderivative(n, a)\n\nReturn the value of the n-th derivative of the polynomial a.\n\n\n\n\n\nderivative(a, r)\n\nPartial differentiation of a::HomogeneousPolynomial series with respect to the r-th variable.\n\n\n\n\n\nderivative(a, r)\n\nPartial differentiation of a::TaylorN series with respect to the r-th variable. The r-th variable may be also specified through its symbol.\n\n\n\n\n\nderivative(a::TaylorN{T}, ntup::NTuple{N,Int})\n\nReturn a TaylorN with the partial derivative of a defined by ntup::NTuple{N,Int}, where the first entry is the number of derivatives with respect to the first variable, the second is the number of derivatives with respect to the second, and so on.\n\n\n\n\n\nderivative(ntup::NTuple{N,Int}, a::TaylorN{T})\n\nReturns the value of the coefficient of a specified by ntup::NTuple{N,Int}, multiplied by the corresponding factorials.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.integrate",
    "page": "API",
    "title": "TaylorSeries.integrate",
    "category": "function",
    "text": "integrate(a, [x])\n\nReturn the integral of a::Taylor1. The constant of integration (0-th order coefficient) is set to x, which is zero if ommitted.\n\n\n\n\n\nintegrate(a, r)\n\nIntegrate the a::HomogeneousPolynomial with respect to the r-th variable. The returned HomogeneousPolynomial has no added constant of integration. If the order of a corresponds to get_order(), a zero HomogeneousPolynomial of 0-th order is returned.\n\n\n\n\n\nintegrate(a, r, [x0])\n\nIntegrate the a::TaylorN series with respect to the r-th variable, where x0 the integration constant and must be independent of the r-th variable; if x0 is ommitted, it is taken as zero.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.gradient",
    "page": "API",
    "title": "TaylorSeries.gradient",
    "category": "function",
    "text": "    gradient(f)\n    ∇(f)\n\nCompute the gradient of the polynomial f::TaylorN.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.jacobian",
    "page": "API",
    "title": "TaylorSeries.jacobian",
    "category": "function",
    "text": "    jacobian(vf)\n    jacobian(vf, [vals])\n\nCompute the jacobian matrix of vf, a vector of TaylorN polynomials, evaluated at the vector vals. If vals is ommited, it is evaluated at zero.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.jacobian!",
    "page": "API",
    "title": "TaylorSeries.jacobian!",
    "category": "function",
    "text": "    jacobian!(jac, vf)\n    jacobian!(jac, vf, [vals])\n\nCompute the jacobian matrix of vf, a vector of TaylorN polynomials evaluated at the vector vals, and write results to jac. If vals is ommited, it is evaluated at zero.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.hessian",
    "page": "API",
    "title": "TaylorSeries.hessian",
    "category": "function",
    "text": "    hessian(f)\n    hessian(f, [vals])\n\nReturn the hessian matrix (jacobian of the gradient) of f::TaylorN, evaluated at the vector vals. If vals is ommited, it is evaluated at zero.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.hessian!",
    "page": "API",
    "title": "TaylorSeries.hessian!",
    "category": "function",
    "text": "    hessian!(hes, f)\n    hessian!(hes, f, [vals])\n\nReturn the hessian matrix (jacobian of the gradient) of f::TaylorN, evaluated at the vector vals, and write results to hes. If vals is ommited, it is evaluated at zero.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.inverse",
    "page": "API",
    "title": "TaylorSeries.inverse",
    "category": "function",
    "text": "inverse(f)\n\nReturn the Taylor expansion of f^-1(t), of order N = f.order, for f::Taylor1 polynomial if the first coefficient of f is zero. Otherwise, an ArgumentError is thrown.\n\nThe algorithm implements Lagrange inversion at t=0 if f(0)=0:\n\nbeginequation*\nf^-1(t) = sum_n=1^N fract^nn left\n    fracrm d^n-1rm d z^n-1left(fraczf(z)right)^n\n    rightvert_z=0\nendequation*\n\n\n\n"
},

{
    "location": "api.html#Base.abs",
    "page": "API",
    "title": "Base.abs",
    "category": "function",
    "text": "abs(a)\n\nReturns a if constant_term(a) > 0 and -a if constant_term(a) < 0 for a <:Union{Taylor1,TaylorN}. Notice that typeof(abs(a)) <: AbstractSeries.\n\n\n\n"
},

{
    "location": "api.html#LinearAlgebra.norm",
    "page": "API",
    "title": "LinearAlgebra.norm",
    "category": "function",
    "text": "norm(x::AbstractSeries, p::Real)\n\nReturns the p-norm of an x::AbstractSeries, defined by\n\nbeginequation*\nleftVert x rightVert_p =  left( sum_k  x_k ^p right)^frac1p\nendequation*\n\nwhich returns a non-negative number.\n\n\n\n"
},

{
    "location": "api.html#Base.isapprox",
    "page": "API",
    "title": "Base.isapprox",
    "category": "function",
    "text": "isapprox(x::AbstractSeries, y::AbstractSeries; rtol::Real=sqrt(eps), atol::Real=0, nans::Bool=false)\n\nInexact equality comparison between polynomials: returns true if norm(x-y,1) <= atol + rtol*max(norm(x,1), norm(y,1)), where x and y are polynomials. For more details, see Base.isapprox.\n\n\n\n\n\n"
},

{
    "location": "api.html#Base.isfinite",
    "page": "API",
    "title": "Base.isfinite",
    "category": "function",
    "text": "isfinite(x::AbstractSeries) -> Bool\n\nTest whether the coefficients of the polynomial x are finite.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.displayBigO",
    "page": "API",
    "title": "TaylorSeries.displayBigO",
    "category": "function",
    "text": "displayBigO(d::Bool) --> nothing\n\nSet/unset displaying of the big 𝒪 notation in  the output of Taylor1 and TaylorN polynomials. The initial value is true.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.use_show_default",
    "page": "API",
    "title": "TaylorSeries.use_show_default",
    "category": "function",
    "text": "use_Base_show(d::Bool) --> nothing\n\nUse Base.show_default method (default show method in Base), or a custom display. The initial value is false, so customized display is used.\n\n\n\n\n\n"
},

{
    "location": "api.html#Functions-and-methods-1",
    "page": "API",
    "title": "Functions and methods",
    "category": "section",
    "text": "Taylor1(::Type{T}, ::Int64=1) where {T<:Number}\nHomogeneousPolynomial(::Type{T}, ::Int) where {T<:Number}\nTaylorN(::Type{T}, ::Int; ::Int=get_order()) where {T<:Number}\nset_variables\nget_variables\nshow_params_TaylorN\nshow_monomials\ngetcoeff\nevaluate\nevaluate!\ntaylor_expand\nupdate!\nderivative\nintegrate\ngradient\njacobian\njacobian!\nhessian\nhessian!\ninverse\nabs\nnorm\nisapprox\nisfinite\ndisplayBigO\nuse_show_default"
},

{
    "location": "api.html#TaylorSeries.ParamsTaylorN",
    "page": "API",
    "title": "TaylorSeries.ParamsTaylorN",
    "category": "type",
    "text": "ParamsTaylorN\n\nDataType holding the current parameters for TaylorN and HomogeneousPolynomial.\n\nFields:\n\norder            :: Int  Order (degree) of the polynomials\nnum_vars         :: Int  Number of variables\nvariable_names   :: Vector{String} Names of the variables\nvariable_symbols :: Vector{Symbol}  Symbols of the variables\n\nThese parameters can be changed using set_variables\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries._InternalMutFuncs",
    "page": "API",
    "title": "TaylorSeries._InternalMutFuncs",
    "category": "type",
    "text": "_InternalMutFuncs\n\nContains parameters and expressions that allow a simple programatic construction for calling the internal mutating functions.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.generate_tables",
    "page": "API",
    "title": "TaylorSeries.generate_tables",
    "category": "function",
    "text": "generate_tables(num_vars, order)\n\nReturn the hash tables coeff_table, index_table, size_table and pos_table. Internally, these are treated as const.\n\nHash tables\n\ncoeff_table :: Array{Array{Array{Int64,1},1},1}\n\nThe i+1-th component contains a vector with the vectors of all the possible combinations of monomials of a HomogeneousPolynomial of order i.\n\nindex_table :: Array{Array{Int64,1},1}\n\nThe i+1-th component contains a vector of (hashed) indices that represent the distinct monomials of a HomogeneousPolynomial of order (degree) i.\n\nsize_table :: Array{Int64,1}\n\nThe i+1-th component contains the number of distinct monomials of the HomogeneousPolynomial of order i, equivalent to length(coeff_table[i]).\n\npos_table :: Array{Dict{Int64,Int64},1}\n\nThe i+1-th component maps the hash index to the (lexicographic) position of the corresponding monomial in coeffs_table.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.generate_index_vectors",
    "page": "API",
    "title": "TaylorSeries.generate_index_vectors",
    "category": "function",
    "text": "generate_index_vectors(num_vars, degree)\n\nReturn a vector of index vectors with num_vars (number of variables) and degree.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.in_base",
    "page": "API",
    "title": "TaylorSeries.in_base",
    "category": "function",
    "text": "in_base(order, v)\n\nConvert vector v of non-negative integers to base order+1.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.make_inverse_dict",
    "page": "API",
    "title": "TaylorSeries.make_inverse_dict",
    "category": "function",
    "text": "make_inverse_dict(v)\n\nReturn a Dict with the enumeration of v: the elements of v point to the corresponding index.\n\nIt is used to construct pos_table from index_table.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.resize_coeffs1!",
    "page": "API",
    "title": "TaylorSeries.resize_coeffs1!",
    "category": "function",
    "text": "resize_coeffs1!{T<Number}(coeffs::Array{T,1}, order::Int)\n\nIf the length of coeffs is smaller than order+1, it resizes coeffs appropriately filling it with zeros.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.resize_coeffsHP!",
    "page": "API",
    "title": "TaylorSeries.resize_coeffsHP!",
    "category": "function",
    "text": "resize_coeffsHP!{T<Number}(coeffs::Array{T,1}, order::Int)\n\nIf the length of coeffs is smaller than the number of coefficients correspondinf to order (given by size_table[order+1]), it resizes coeffs appropriately filling it with zeros.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.constant_term",
    "page": "API",
    "title": "TaylorSeries.constant_term",
    "category": "function",
    "text": "constant_term(a)\n\nReturn the constant value (zero order coefficient) for Taylor1 and TaylorN.\n\n\n\n\n\n"
},

{
    "location": "api.html#LinearAlgebra.mul!",
    "page": "API",
    "title": "LinearAlgebra.mul!",
    "category": "function",
    "text": "mul!(c, a, b, k::Int) --> nothing\n\nUpdate the k-th expansion coefficient c[k] of c = a * b, where all c, a, and b are either Taylor1 or TaylorN.\n\nThe coefficients are given by\n\nc_k = sum_j=0^k a_j b_k-j\n\n\n\n"
},

{
    "location": "api.html#LinearAlgebra.mul!-Tuple{HomogeneousPolynomial,HomogeneousPolynomial,HomogeneousPolynomial}",
    "page": "API",
    "title": "LinearAlgebra.mul!",
    "category": "method",
    "text": "mul!(c, a, b) --> nothing\n\nReturn c = a*b with no allocation; all arguments are HomogeneousPolynomial.\n\n\n\n\n\n"
},

{
    "location": "api.html#LinearAlgebra.mul!-Union{Tuple{T}, Tuple{Array{Taylor1{T},1},Union{Array{T,2}, SparseMatrixCSC{T,Ti} where Ti<:Integer},Array{Taylor1{T},1}}} where T<:Number",
    "page": "API",
    "title": "LinearAlgebra.mul!",
    "category": "method",
    "text": "mul!(Y, A, B)\n\nMultiply A*B and save the result in Y.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.div!",
    "page": "API",
    "title": "TaylorSeries.div!",
    "category": "function",
    "text": "div!(c, a, b, k::Int, ordfact::Int=0)\n\nCompute the k-th expansion coefficient c[k] of c = a / b, where all c, a and b are either Taylor1 or TaylorN.\n\nThe coefficients are given by\n\nc_k =  frac1b_0 big(a_k - sum_j=0^k-1 c_j b_k-jbig)\n\nFor Taylor1 polynomials, ordfact is the order of the factorized term of the denominator.\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.pow!",
    "page": "API",
    "title": "TaylorSeries.pow!",
    "category": "function",
    "text": "pow!(c, a, r::Real, k::Int, k0::Int=0)\n\nUpdate the k-th expansion coefficient c[k] of c = a^r, for both c and a either Taylor1 or TaylorN.\n\nThe coefficients are given by\n\nc_k = frac1k a_0 sum_j=0^k-1 big(r(k-j) -jbig)a_k-j c_j\n\nFor Taylor1 polynomials, k0 is the order of the first non-zero coefficient of a.\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.square",
    "page": "API",
    "title": "TaylorSeries.square",
    "category": "function",
    "text": "square(a::AbstractSeries) --> typeof(a)\n\nReturn a^2; see TaylorSeries.sqr!.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.sqr!",
    "page": "API",
    "title": "TaylorSeries.sqr!",
    "category": "function",
    "text": "sqr!(c, a, k::Int) --> nothing\n\nUpdate the k-th expansion coefficient c[k] of c = a^2, for both c and a either Taylor1 or TaylorN.\n\nThe coefficients are given by\n\nbegineqnarray*\nc_k  =  2 sum_j=0^(k-1)2 a_k-j a_j\n    text if k is odd \nc_k  =  2 sum_j=0^(k-2)2 a_k-j a_j + (a_k2)^2\n    text if k is even \nendeqnarray*\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.sqr!-Union{Tuple{T}, Tuple{HomogeneousPolynomial{T},HomogeneousPolynomial{T}}} where T<:Union{Real, Complex, Taylor1}",
    "page": "API",
    "title": "TaylorSeries.sqr!",
    "category": "method",
    "text": "sqr!(c, a)\n\nReturn c = a*a with no allocation; all parameters are HomogeneousPolynomial.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.sqrt!",
    "page": "API",
    "title": "TaylorSeries.sqrt!",
    "category": "function",
    "text": "sqrt!(c, a, k::Int, k0::Int=0)\n\nCompute the k-th expansion coefficient c[k] of c = sqrt(a) for bothc and a either Taylor1 or TaylorN.\n\nThe coefficients are given by\n\nbegineqnarray*\nc_k = frac12 c_0 big( a_k - 2 sum_j=1^(k-1)2 c_k-jc_jbig)\n    text if k is odd \nc_k = frac12 c_0 big( a_k - 2 sum_j=1^(k-2)2 c_k-jc_j\n    - (c_k2)^2big) text if k is even\nendeqnarray*\n\nFor Taylor1 polynomials, k0 is the order of the first non-zero coefficient, which must be even.\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.exp!",
    "page": "API",
    "title": "TaylorSeries.exp!",
    "category": "function",
    "text": "exp!(c, a, k) --> nothing\n\nUpdate the k-th expansion coefficient c[k+1] of c = exp(a) for both c and a either Taylor1 or TaylorN.\n\nThe coefficients are given by\n\nbeginequation*\nc_k = frac1k sum_j=0^k-1 (k-j) a_k-j c_j\nendequation*\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.log!",
    "page": "API",
    "title": "TaylorSeries.log!",
    "category": "function",
    "text": "log!(c, a, k) --> nothing\n\nUpdate the k-th expansion coefficient c[k+1] of c = log(a) for both c and a either Taylor1 or TaylorN.\n\nThe coefficients are given by\n\nbeginequation*\nc_k = frac1a_0 big(a_k - frac1k sum_j=0^k-1 j a_k-j c_j big)\nendequation*\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.sincos!",
    "page": "API",
    "title": "TaylorSeries.sincos!",
    "category": "function",
    "text": "sincos!(s, c, a, k) --> nothing\n\nUpdate the k-th expansion coefficients s[k+1] and c[k+1] of s = sin(a) and c = cos(a) simultaneously, for s, c and a either Taylor1 or TaylorN.\n\nThe coefficients are given by\n\nbegineqnarray*\ns_k =  frac1ksum_j=0^k-1 (k-j) a_k-j c_j \nc_k = -frac1ksum_j=0^k-1 (k-j) a_k-j s_j\nendeqnarray*\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.tan!",
    "page": "API",
    "title": "TaylorSeries.tan!",
    "category": "function",
    "text": "tan!(c, a, p, k::Int) --> nothing\n\nUpdate the k-th expansion coefficients c[k+1] of c = tan(a), for c and a either Taylor1 or TaylorN; p = c^2 and is passed as an argument for efficiency.\n\nThe coefficients are given by\n\nbeginequation*\nc_k = a_k + frac1k sum_j=0^k-1 (k-j) a_k-j p_j\nendequation*\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.asin!",
    "page": "API",
    "title": "TaylorSeries.asin!",
    "category": "function",
    "text": "asin!(c, a, r, k)\n\nUpdate the k-th expansion coefficients c[k+1] of c = asin(a), for c and a either Taylor1 or TaylorN; r = sqrt(1-c^2) and is passed as an argument for efficiency.\n\nbeginequation*\nc_k = frac1 sqrtr_0 \n    big( a_k - frac1k sum_j=1^k-1 j r_k-j c_j big)\nendequation*\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.acos!",
    "page": "API",
    "title": "TaylorSeries.acos!",
    "category": "function",
    "text": "acos!(c, a, r, k)\n\nUpdate the k-th expansion coefficients c[k+1] of c = acos(a), for c and a either Taylor1 or TaylorN; r = sqrt(1-c^2) and is passed as an argument for efficiency.\n\nbeginequation*\nc_k = - frac1 r_0 \n    big( a_k - frac1k sum_j=1^k-1 j r_k-j c_j big)\nendequation*\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.atan!",
    "page": "API",
    "title": "TaylorSeries.atan!",
    "category": "function",
    "text": "atan!(c, a, r, k)\n\nUpdate the k-th expansion coefficients c[k+1] of c = atan(a), for c and a either Taylor1 or TaylorN; r = 1+a^2 and is passed as an argument for efficiency.\n\nbeginequation*\nc_k = frac1r_0big(a_k - frac1k sum_j=1^k-1 j r_k-j c_jbig)\nendequation*\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.sinhcosh!",
    "page": "API",
    "title": "TaylorSeries.sinhcosh!",
    "category": "function",
    "text": "sinhcosh!(s, c, a, k)\n\nUpdate the k-th expansion coefficients s[k+1] and c[k+1] of s = sinh(a) and c = cosh(a) simultaneously, for s, c and a either Taylor1 or TaylorN.\n\nThe coefficients are given by\n\nbegineqnarray*\ns_k = frac1k sum_j=0^k-1 (k-j) a_k-j c_j \nc_k = frac1k sum_j=0^k-1 (k-j) a_k-j s_j\nendeqnarray*\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.tanh!",
    "page": "API",
    "title": "TaylorSeries.tanh!",
    "category": "function",
    "text": "tanh!(c, a, p, k)\n\nUpdate the k-th expansion coefficients c[k+1] of c = tanh(a), for c and a either Taylor1 or TaylorN; p = a^2 and is passed as an argument for efficiency.\n\nbeginequation*\nc_k = a_k - frac1k sum_j=0^k-1 (k-j) a_k-j p_j\nendequation*\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries.derivative!",
    "page": "API",
    "title": "TaylorSeries.derivative!",
    "category": "function",
    "text": "derivative!(res, a) --> nothing\n\nIn-place version of derivative. Compute the Taylor1 polynomial of the differential of a::Taylor1 and save it into res. The last coefficient is set to zero.\n\n\n\n\n\nderivative!(p, a, k) --> nothing\n\nUpdate in-place the k-th expansion coefficient p[k] of p = derivative(a) for both p and a Taylor1.\n\nThe coefficients are given by\n\np_k = (k+1)a_k+1\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries._internalmutfunc_call",
    "page": "API",
    "title": "TaylorSeries._internalmutfunc_call",
    "category": "function",
    "text": "_internalmutfunc_call( fn :: _InternalMutFuncs )\n\nCreates the appropriate call to the internal mutating function defined by the _InternalMutFuncs object. This is used to construct _dict_unary_calls and _dict_binary_calls. The call contains the prefix TaylorSeries..\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries._dict_unary_ops",
    "page": "API",
    "title": "TaylorSeries._dict_unary_ops",
    "category": "constant",
    "text": "_dict_binary_ops\n\nDict{Symbol, Array{Any,1}} with the information to construct the _InternalMutFuncs related to unary operations.\n\nThe keys correspond to the function symbols.\n\nThe arguments of the array are the function name (e.g. add!), a tuple with the function arguments, and an Expr with the calling pattern. The convention for the arguments of the functions and the calling pattern is to use :_res for the (mutated) result, :_arg1, for the required argument, possibly :_aux when there is an auxiliary expression needed, and :_k for the computed order of :_res. When an auxiliary expression is required, and Expr defining its calling pattern is added as the last entry of the vector.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries._dict_binary_calls",
    "page": "API",
    "title": "TaylorSeries._dict_binary_calls",
    "category": "constant",
    "text": "_dict_binary_calls::Dict{Symbol, NTuple{2,Expr}}\n\nDictionary with the expressions that define the internal binary functions and the auxiliary functions, whenever they exist. The keys correspond to those functions, passed as symbols, with the defined internal mutating functions.\n\nEvaluating the entries generates symbols that represent the actual calls to the internal mutating functions.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries._dict_unary_calls",
    "page": "API",
    "title": "TaylorSeries._dict_unary_calls",
    "category": "constant",
    "text": "_dict_unary_calls::Dict{Symbol, NTuple{2,Expr}}\n\nDictionary with the expressions that define the internal unary functions and the auxiliary functions, whenever they exist. The keys correspond to those functions, passed as symbols, with the defined internal mutating functions.\n\nEvaluating the entries generates expressions that represent the actual calls to the internal mutating functions.\n\n\n\n\n\n"
},

{
    "location": "api.html#TaylorSeries._dict_binary_ops",
    "page": "API",
    "title": "TaylorSeries._dict_binary_ops",
    "category": "constant",
    "text": "_dict_binary_ops\n\nDict{Symbol, Array{Any,1}} with the information to construct the _InternalMutFuncs related to binary operations.\n\nThe keys correspond to the function symbols.\n\nThe arguments of the array are the function name (e.g. add!), a tuple with the function arguments, and an Expr with the calling pattern. The convention for the arguments of the functions and the calling pattern is to use :_res for the (mutated) result, :_arg1 and _arg2 for the required arguments, and :_k for the computed order of :_res.\n\n\n\n\n\n"
},

{
    "location": "api.html#Internals-1",
    "page": "API",
    "title": "Internals",
    "category": "section",
    "text": "ParamsTaylorN\n_InternalMutFuncs\ngenerate_tables\ngenerate_index_vectors\nin_base\nmake_inverse_dict\nresize_coeffs1!\nresize_coeffsHP!\nconstant_term\nmul!\nmul!(::HomogeneousPolynomial, ::HomogeneousPolynomial, ::HomogeneousPolynomial)\nmul!(::Vector{Taylor1{T}}, ::Union{Matrix{T},SparseMatrixCSC{T}},::Vector{Taylor1{T}}) where {T<:Number}\ndiv!\npow!\nsquare\nsqr!\nsqr!(::HomogeneousPolynomial{T}, ::HomogeneousPolynomial{T}) where {T<:NumberNotSeriesN}\nsqrt!\nexp!\nlog!\nsincos!\ntan!\nasin!\nacos!\natan!\nsinhcosh!\ntanh!\nderivative!\n_internalmutfunc_call\n_dict_unary_ops\n_dict_binary_calls\n_dict_unary_calls\n_dict_binary_ops"
},

{
    "location": "api.html#Index-1",
    "page": "API",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"api.md\"]\nModule = [\"TaylorSeries\"]\nOrder = [:type, :function]"
},

]}
