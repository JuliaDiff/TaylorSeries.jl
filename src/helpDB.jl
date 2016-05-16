# This file is part of the TaylorSeries.jl Julia package, MIT license
#
#`docstrings` for the Taylor1.jl and TaylorN.jl scripts
#


                                        ### Exclusive Operations ###



## Constructors ##
@doc """
    Taylor1{T<:Number} <: Number

DataType for polynomial expansions in one independent variable.

**Fields:**

- `coeffs :: Array{T,1}` Expansion coefficients; the \$i\$-th
component is the coefficient of degree \$i-1\$ of the expansion.
- `order  :: Int64` Maximum order (degree) of the polynomial.
"""->
Taylor1

@doc """
    taylor1_variable(T, [order=1])
    taylor1_variable([order=1])

Short-cut to define the independent variable as a `Taylor1` polynomial of
given `order`. If `T::Type` is ommitted, `Float64` is assumend.
"""-> 
taylor1_variable

## get_coeff ##
@doc """
    get_coeff(a, n)

Return the coefficient of order `n::Int` of a `a::Taylor1` polynomial.
""" ->
get_coeff

## Differentiating ##
doc"""
    diffTaylor(a)

Return the `Taylor1` polynomial of the differential of `a::Taylor1`. 

The last coefficient is set to zero.
"""
diffTaylor

## Integrating ##
doc"""
    integTaylor(a, x)
    integTaylor(a)

Returns the integral of `a::Taylor1`. 

The constant of integration (0th order coefficient) is set to `x`, which is zero if ommitted.
"""
integTaylor

## Evaluating ##
doc"""
    evaluate(a, dx)
    evaluate(a)

Evaluates a `Taylor1` polynomial using Horner's rule (hand coded).


    evaluate(a, x)

Returns the substitution of `x::Taylor1` as independent variable in
`a::Taylor1`.
"""
evaluate

## Returns de n-th derivative of a series expansion
doc"""
    deriv(a, [n=1])

Returns the value of the `n`-th derivative of `a`.
"""
deriv


                                        ### Function Extensions ###


## Exp ##
doc"""
    exp(a)

Compute $e^a$'s expansion of order `a.order` as a `Taylor1` object for `a::Taylor1`.
"""
exp

## Sin ##
doc"""
    sin(a)

Compute $\sin{a}$'s expansion of order `a.order` as a `Taylor1` object for `a::Taylor1`.
"""
sin

## Cos ## 
doc"""
    cos(a)

Computes $\cos{a}$'s expansion of order `a.order` as a `Taylor1` object for `a::Taylor1`.
"""
cos

## Log ##
doc"""
    log(a)

Computes the $\log{a}$'s expansion of order `a.order` as a `Taylor1` object for `a:Taylor1`.
"""
log

## Tan ##
doc"""
    tan(a)

Computes the $\tan{a}$'s expansion of order `a.order` as a `Taylor1` object for `a::Taylor1`.
"""
tan

## abs function ##
doc"""
    abs(a)

Computes `a` or `-a` depending on the 0-th order coefficient
of `a::Taylor1`.

<br>
If it is zero, an `ArgumentError` is thrown.
"""
abs

## power ##
doc"""
    ^(a, x)

Computes $a^x$ as a `Taylor1` object with `a::Taylor1` and `x::Number`.

<br>
If `x::Real` and the 0th order coefficient is non-zero, an `ArgumentError` is thrown.
"""
^

## Square root ##
doc"""
    sqrt(a)

Computes $\sqrt{a}$'s expansion of order `a.order` as an `Taylor1` object.

If the first non-vanishing coefficient of `a` is an **odd power**, and `ArgumentError` will be thrown.  
"""
sqrt


                                        ### Operation Coefficients ###



# Homogeneous coefficient for the multiplication
doc"""
Computes the **`k-th` expansion coefficient** of $p(x) = f(x) g(x)$ as

<center>
$ p_k = \sum_{j=0}^k f_j g_{k-j} $
</center>

with $f(x)$  and $g(x)$ analitical functions.
"""
mulHomogCoef

# Homogeneous coefficient for the division
doc"""
Computes the **`k-th` expansion coefficient** of $d(x) = f(x) / g(x)$ as 

<center>
$ d_k =  \frac{1}{g_0} ( f_r - \sum_{j=0}^{k-1} d_j g_{k-j} ) $
</center>

with $f(x)$  and $g(x)$ analitical functions.
"""
divHomogCoef

# Homogeneous coefficients for real power
doc"""
Computes the **`k-th` expansion coefficient** of $p(x) = f(x)^β$ as 

<center>
$ p_k = \frac{1}{k f_0} \sum_{j=0}^{k-1} (β(k-j) -j) - j)f_{r-j} p_j $
</center>

with $f(x)$ an analitical function.
"""
powHomogCoef

# Homogeneous coefficients for square
doc"""
Computes the **`k-th` expansion coefficient** of $s(x) = f(x)^2$ as

<center>
$ s_k = 2 \sum_{j=0}^{(k-1)/2} f_{k-j} f_j $
</center>

when `k`is **odd**,

<center>
$ s_k = 2 \sum_{j=0}^{(k-2)/2} ( f_{k-j} f_j + (f_{r/2})^2 )$ 
</center>

when `r`is **even**,

with $f(x)$ an analitical function.
"""
squareHomogCoef

# Homogeneous coefficients for the square-root
doc"""
Computes the **`k-th` expansion coefficient** of $s(x) = \sqrt{f(x)}$ as

<center>
$ s_k = \frac{1}{2 s_0} ( f_k - 2 \sum_{j=0}^{(k-1)/2} s_{k-j}s_j ) $
</center>

when `k` is **odd**,

<center>
$ s_k = \frac{1}{2 s_0} ( f_k - 2 \sum_{j=0}^{(k2-)/2} ( s_{k-j}s_j - (s_{k/2})^2 ) $
</center>

when `k`is **even**,

with $f(x)$ an analitical function.
"""
sqrtHomogCoef

# Homogeneous coefficients for exp
doc"""
Computes the **`k-th` expansion coefficient** of $e(x) = \exp{f(x)}$ as

<center>
$ e_k = \frac{1}{k} \sum_{j=0}^{k-1} (k-j) f_{k-j} e_j $
</center>

with $f(x)$ an analitical function.
"""
expHomogCoef

# Homogeneous coefficients for log
doc"""
Computes the **`k-th` expansion coefficient** of $l(x) = \log{f(x)}$ as

<center>
$ l_k = \frac{1}{f_0} ( f_k - \frac{1}{k} \sum_{j=0}^{k-1} j f_{k-j} l_j ) $
</center>

with $f(x)$ an analitical function.
"""
logHomogCoef

# Homogeneous coefficients for sincos
doc"""
Computes the **`k-th` expansion coefficients ** of $s(x) = \sin{f(x)}$ and $c(x) = \cos{f(x)}$ as

<center>
$ s_k = \frac{1}{k} \sum_{j=0}^{k-1} (k-j) f_{k-j} c_j $
</center>

and

<center>
$ c_k = -\frac{1}{k} \sum_{j=0}^{k-1} (k-j) f_{k-j} s_j $
</center>

with $f(x)$ an analitical function.
"""
sincosHomogCoef

# Homogeneous coefficients for tan
doc"""
Computes the **`r-th` expansion coefficient** of $t(x) = \tan{f(x)}$ as 

<center>
$ t_k = f_k + \frac{1}{k} \sum_{j=0}^{k-1} (k-j) f_{k-j} p_j $
</center>

with $f(x)$ an analitical function and $p(x) := t(x)^2$.
"""
tanHomogCoef