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
JetSpace
```

## Functions and methods

```@docs
Taylor1(::Type{T}, ::Int) where {T<:Number}
HomogeneousPolynomial(::Type{T}, ::Int) where {T<:Number}
TaylorN(::Type{T}, ::Int; ::Int=get_order()) where {T<:Number}
set_variables
variables
get_order
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
jacobianmatrix
hessian
hessian!
hessianmatrix
constant_term
linear_polynomial
nonlinear_polynomial
inverse
inverse_map
abs
norm
isapprox
isless
isfinite
displayBigO
use_show_default
set_taylor1_varname
space
```

## Internals

```@docs
ParamsTaylor1
ParamsTaylorN
HomogeneousProductTable
_InternalMutFuncs
generate_tables
_homogeneous_product_table
_init_output_major_product_table!
generate_multiplication_tables
_product_table
_init_product_table!
generate_index_vectors
in_base
make_inverse_dict
_sync_legacy_tables!
set_default_space!
_coeffsHP
_coeffsTN
_check_same_space
numtype
mul!
mul!(::HomogeneousPolynomial, ::HomogeneousPolynomial, ::HomogeneousPolynomial)
mul_scalar!(::HomogeneousPolynomial, ::NumberNotSeries, ::HomogeneousPolynomial, ::HomogeneousPolynomial)
mul!(::Vector{Taylor1{T}}, ::Union{Matrix{T},SparseMatrixCSC{T}},::Vector{Taylor1{T}}) where {T<:Number}
div!
pow!
square
sqr!
accsqr!
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
_isthinzero
_internalmutfunc_call
_dict_unary_ops
_dict_binary_calls
_dict_unary_calls
_dict_binary_ops
_populate_dicts!
@isonethread
```

## Index

```@index
Pages = ["api.md"]
Order = [:type, :function]
```
