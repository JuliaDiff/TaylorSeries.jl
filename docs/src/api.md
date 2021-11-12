# Library

---

```@meta
CurrentModule = TaylorSeries
```

## Module
```@docs
TaylorSeries
```

## Types

```@docs
Taylor1
HomogeneousPolynomial
TaylorN
AbstractSeries
```

## Functions and methods

```@docs
Taylor1(::Type{T}, ::Int) where {T<:Number}
HomogeneousPolynomial(::Type{T}, ::Int) where {T<:Number}
TaylorN(::Type{T}, ::Int; ::Int=get_order()) where {T<:Number}
set_variables
get_variables
show_params_TaylorN
show_monomials
getcoeff
evaluate
evaluate!
taylor_expand
update!
differentiate
derivative
integrate
gradient
jacobian
jacobian!
hessian
hessian!
constant_term
linear_polynomial
inverse
abs
norm
isapprox
isfinite
displayBigO
use_show_default
set_taylor1_varname
```

## Internals

```@docs
ParamsTaylor1
ParamsTaylorN
_InternalMutFuncs
generate_tables
generate_index_vectors
in_base
make_inverse_dict
resize_coeffs1!
resize_coeffsHP!
numtype
mul!
mul!(::HomogeneousPolynomial, ::HomogeneousPolynomial, ::HomogeneousPolynomial)
mul!(::Vector{Taylor1{T}}, ::Union{Matrix{T},SparseMatrixCSC{T}},::Vector{Taylor1{T}}) where {T<:Number}
div!
pow!
square
sqr!
sqr!(::HomogeneousPolynomial{T}, ::HomogeneousPolynomial{T}) where {T<:NumberNotSeriesN}
sqrt!
exp!
log!
sincos!
tan!
asin!
acos!
atan!
sinhcosh!
tanh!
asinh!
acosh!
atanh!
differentiate!
_internalmutfunc_call
_dict_unary_ops
_dict_binary_calls
_dict_unary_calls
_dict_binary_ops
```

## Index

```@index
Pages = ["api.md"]
Module = ["TaylorSeries"]
Order = [:type, :function]
```
