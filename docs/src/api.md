# Library

---

```@meta
CurrentModule = TaylorSeries
```

## Types

```@docs
Taylor1
HomogeneousPolynomial
TaylorN
AbstractSeries
ParamsTaylorN
```

## Functions and methods

```@docs
Taylor1{T<:Number}(::Type{T}, ::Int)
HomogeneousPolynomial{T<:Number}(::Type{T}, ::Int)
TaylorN{T<:Number}(::Type{T}, nv::Int; [order::Int=get_order()])
set_variables
show_params_TaylorN
get_coeff
evaluate
evaluate!
derivative
integrate
gradient
jacobian
jacobian!
hessian
hessian!
inverse
```

## Internals

```@docs
generate_tables
generate_index_vectors
in_base
make_inverse_dict
order_posTb
max_order
resize_coeffs1!
resize_coeffsHP!
zero_korder
constant_term
mul!
mul!(c::HomogeneousPolynomial, a::HomogeneousPolynomial, b::HomogeneousPolynomial)
div!
pow!
square
sqr!
sqr!(c::HomogeneousPolynomial, a::HomogeneousPolynomial)
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
A_mul_B!
```

## Index

```@index
Pages = ["api.md"]
Module = ["TaylorSeries"]
Order = []
```
