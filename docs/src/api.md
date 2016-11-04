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
*
/
^
sqrt
exp
log
sin
cos
tan
asin
acos
atan
sinh
cosh
tanh
gaussian
erf
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
asinHomogCoef
acosHomogCoef
atanHomogCoef
sinhcoshHomogCoef
tanhHomogCoef
gaussHomogCoef
erfHomogCoef
A_mul_B!
```
