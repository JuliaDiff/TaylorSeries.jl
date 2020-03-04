# Examples

---

```@meta
CurrentModule = TaylorSeries
```

## Four-square identity

The first example shows that the four-square identity holds:
```math
(a_1+a_2+a_3+a_4)\cdot(b_1+b_2+b_3+b_4) = \\
  \qquad (a_1 b_1 - a_2 b_2 - a_3 b_3 -a_4 b_4)^2 +
  (a_1 b_2 - a_2 b_1 - a_3 b_4 -a_4 b_3)^2 + \\
  \qquad (a_1 b_3 - a_2 b_4 - a_3 b_1 -a_4 b_2)^2 +
  (a_1 b_4 - a_2 b_3 - a_3 b_2 -a_4 b_1)^2,
```
which was originally proved by Euler. The code can also be found in
[this test](https://github.com/JuliaDiff/TaylorSeries.jl/blob/master/test/identities_Euler.jl) of the package.

First, we reset the maximum degree of the polynomial to 4, since the RHS
of the equation has *a priori* terms of fourth order, and define the 8
independent variables.

```@repl euler
using TaylorSeries
# Define the variables α₁, ..., α₄, β₁, ..., β₄
make_variable(name, index::Int) = string(name, TaylorSeries.subscriptify(index))
variable_names = [make_variable("α", i) for i in 1:4]
append!(variable_names, [make_variable("β", i) for i in 1:4])
# Create the TaylorN variables (order=4, numvars=8)
a1, a2, a3, a4, b1, b2, b3, b4 = set_variables(variable_names, order=4)
a1 # variable a1
```

Now we compute each term appearing in Eq. (\ref{eq:Euler})

```@repl euler
# left-hand side
lhs1 = a1^2 + a2^2 + a3^2 + a4^2 ;
lhs2 = b1^2 + b2^2 + b3^2 + b4^2 ;
lhs = lhs1 * lhs2;
# right-hand side
rhs1 = (a1*b1 - a2*b2 - a3*b3 - a4*b4)^2 ;
rhs2 = (a1*b2 + a2*b1 + a3*b4 - a4*b3)^2 ;
rhs3 = (a1*b3 - a2*b4 + a3*b1 + a4*b2)^2 ;
rhs4 = (a1*b4 + a2*b3 - a3*b2 + a4*b1)^2 ;
rhs = rhs1 + rhs2 + rhs3 + rhs4;
```

We now compare the two sides of the identity,

```@repl euler
lhs == rhs
```

The identity is satisfied. ``\square``



## Fateman test

Richard J. Fateman, from Berkley, proposed as a stringent test
of polynomial multiplication
the evaluation of ``s\cdot(s+1)``, where ``s = (1+x+y+z+w)^{20}``. This is
implemented in
the function `fateman1` below. We shall also consider the form
``s^2+s`` in `fateman2`,
which involves fewer operations (and makes a fairer comparison to what
Mathematica does).

```@repl fateman
using TaylorSeries
const order = 20
const x, y, z, w = set_variables(Int128, "x", numvars=4, order=2order)
function fateman1(degree::Int)
    T = Int128
    s = one(T) + x + y + z + w
    s = s^degree
    s * ( s + one(T) )
end
```

(In the following lines, which are run when the documentation is built,
by some reason the timing appears before the command executed.)

```@repl fateman
@time fateman1(0);
# hide
@time f1 = fateman1(20);
```

Another implementation of the same, but exploiting optimizations
related to `^2` yields:

```@repl fateman
function fateman2(degree::Int)
    T = Int128
    s = one(T) + x + y + z + w
    s = s^degree
    s^2 + s
end
fateman2(0);
@time f2 = fateman2(20); # the timing appears above
```

We note that the above functions use expansions in `Int128`. This is actually
required, since some coefficients are larger than `typemax(Int)`:

```@repl fateman
getcoeff(f2, (1,6,7,20)) # coefficient of x y^6 z^7 w^{20}
ans > typemax(Int)
length(f2)
sum(TaylorSeries.size_table)
set_variables("x", numvars=2, order=6) # hide
```

These examples show that
`fateman2` is nearly twice as fast as `fateman1`, and that the series has 135751
monomials in 4 variables.


### Bechmarks

The functions described above have been compared against Mathematica v11.1.
The relevant files used for benchmarking can be found
[here](https://github.com/JuliaDiff/TaylorSeries.jl/tree/master/perf).
Running on a MacPro with Intel-Xeon processors 2.7GHz, we obtain that
Mathematica requires on average (5 runs) 3.075957 seconds for the computation,
while for `fateman1` and `fateman2` above we obtain 2.15408 and 1.08337,
respectively.

Then, with the current version of `TaylorSeries.jl` and using Julia v0.7.0,
our implementation of `fateman1` is about 30%-40% faster.
(The original test by Fateman corresponds to `fateman1` above, which
avoids some optimizations related to squaring; the implementation in Mathematica
is done such that this optimization does not occur.)
