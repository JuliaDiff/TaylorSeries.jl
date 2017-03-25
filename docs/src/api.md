# Library

---

```@meta
CurrentModule = TaylorSeries
```


### Types

```@docs
Taylor1
HomogeneousPolynomial
TaylorN
ParamsTaylorN
```

### Functions and methods

```@docs
Taylor1{T<:Number}(::Type{T},::Int)
TaylorN{T<:Number}(::Type{T}, nv::Int)
set_variables
show_params_TaylorN
get_coeff
evaluate
evaluate!
derivative
integrate
gradient
jacobian
hessian
```

### Internals

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
div!
pow!
sqr!
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
