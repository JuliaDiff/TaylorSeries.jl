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
taylorN_variable
set_variables
show_params_TaylorN
get_coeff
evaluate
derivative
integrate
gradient
jacobian
hessian
*
/
^
sqrt
exp
log
sin
cos
tan
abs
```

### Internals

```@docs
generate_tables
generate_index_vectors
in_base
make_inverse_dict
order_posTb
mul!
mulHomogCoef
divHomogCoef
powHomogCoef
squareHomogCoef
sqrtHomogCoef
expHomogCoef
logHomogCoef
sincosHomogCoef
tanHomogCoef
A_mul_B!
```
