# Examples

## 1. Four-square identity

The first example shows that the four-square identity holds:
\\begin{eqnarray}
(a_1+a_2+a_3+a_4)\\cdot(b_1+b_2+b_3+b_4) & = &
     (a_1 b_1 - a_2 b_2 - a_3 b_3 -a_4 b_4)^2 + \\qquad \\nonumber \\\\
\\label{eq:Euler}
  & & (a_1 b_2 - a_2 b_1 - a_3 b_4 -a_4 b_3)^2 + \\\\
  & & (a_1 b_3 - a_2 b_4 - a_3 b_1 -a_4 b_2)^2 + \\nonumber \\\\
  & & (a_1 b_4 - a_2 b_3 - a_3 b_2 -a_4 b_1)^2, \\nonumber
\\end{eqnarray}
as proved by Euler. The code can we found in one of the tests of the package.

First, we reset the maximum degree of the polynomial to 4, since the RHS
of the equation
has *a priori* terms of fourth order, and the number of independent variables to
8.

    {meta}
    DocTestSetup = quote
        using TaylorSeries
        a1, a2, a3, a4, b1, b2, b3, b4 = set_variables("α₁ α₂ α₃ α₄ β₁ β₂ β₃ β₄", order=4)
    end

```julia
julia> # Define the variables α₁, ..., α₄, β₁, ..., β₄
       make_variable(name, index::Int) = string(name, TaylorSeries.subscriptify(index))
make_variable (generic function with 1 method)

julia> variable_names = [make_variable("α", i) for i in 1:4];

julia> append!(variable_names, [make_variable("β", i) for i in 1:4]);

julia> # Create the Taylor objects (order 4, numvars=8)
       a1, a2, a3, a4, b1, b2, b3, b4 = set_variables(variable_names, order=4);

julia> a1
 1.0 α₁ + 𝒪(‖x‖⁵)

julia> b1
 1.0 β₁ + 𝒪(‖x‖⁵)

```

Now we compute each term appearing in (\\ref{eq:Euler}), and compare them

```julia
julia> # left-hand side

julia> lhs1 = a1^2 + a2^2 + a3^2 + a4^2 ;

julia> lhs2 = b1^2 + b2^2 + b3^2 + b4^2 ;

julia> lhs = lhs1 * lhs2
 1.0 α₁² β₁² + 1.0 α₁² β₂² + 1.0 α₁² β₃² + 1.0 α₁² β₄² + 1.0 α₂² β₁² + 1.0 α₂² β₂² + 1.0 α₂² β₃² + 1.0 α₂² β₄² + 1.0 α₃² β₁² + 1.0 α₃² β₂² + 1.0 α₃² β₃² + 1.0 α₃² β₄² + 1.0 α₄² β₁² + 1.0 α₄² β₂² + 1.0 α₄² β₃² + 1.0 α₄² β₄² + 𝒪(‖x‖⁵)

julia> # right-hand side

julia> rhs1 = (a1*b1 - a2*b2 - a3*b3 - a4*b4)^2 ;

julia> rhs2 = (a1*b2 + a2*b1 + a3*b4 - a4*b3)^2 ;

julia> rhs3 = (a1*b3 - a2*b4 + a3*b1 + a4*b2)^2 ;

julia> rhs4 = (a1*b4 + a2*b3 - a3*b2 + a4*b1)^2
 1.0 α₁² β₄² + 2.0 α₁ α₂ β₃ β₄ - 2.0 α₁ α₃ β₂ β₄ + 2.0 α₁ α₄ β₁ β₄ + 1.0 α₂² β₃² - 2.0 α₂ α₃ β₂ β₃ + 2.0 α₂ α₄ β₁ β₃ + 1.0 α₃² β₂² - 2.0 α₃ α₄ β₁ β₂ + 1.0 α₄² β₁² + 𝒪(‖x‖⁵)

julia> rhs = rhs1 + rhs2 + rhs3 + rhs4
 1.0 α₁² β₁² + 1.0 α₁² β₂² + 1.0 α₁² β₃² + 1.0 α₁² β₄² + 1.0 α₂² β₁² + 1.0 α₂² β₂² + 1.0 α₂² β₃² + 1.0 α₂² β₄² + 1.0 α₃² β₁² + 1.0 α₃² β₂² + 1.0 α₃² β₃² + 1.0 α₃² β₄² + 1.0 α₄² β₁² + 1.0 α₄² β₂² + 1.0 α₄² β₃² + 1.0 α₄² β₄² + 𝒪(‖x‖⁵)

julia> lhs == rhs
true

```

The identity is thus satisfied. $\\square$.



## 2. Fateman test

Richard J. Fateman, from Berkley, proposed as a stringent test
of polynomial multiplication
the evaluation of $s*(s+1)$, where $s = (1+x+y+z+w)^{20}$. This is
implemented in
the function `fateman1`. We shall also evaluate the form $s^2+s$ in `fateman2`,
which involves fewer operations (and makes a fairer comparison to what
Mathematica does). Below we have used Julia v0.4.

    {meta}
    DocTestSetup = quote
        using TaylorSeries
        set_variables("x", numvars=4, order=40);
    end

```julia
julia> set_variables("x", numvars=4, order=40);

julia> function fateman1(degree::Int)
           T = Int128
           oneH = HomogeneousPolynomial(one(T), 0)
           # s = 1 + x + y + z + w
           s = TaylorN( [oneH, HomogeneousPolynomial([one(T),one(T),one(T),one(T)],1)], degree )
           s = s^degree
           # s is converted to order 2*ndeg
           s = TaylorN(s, 2*degree)
           s * ( s+TaylorN(oneH, 2*degree) )
       end
fateman1 (generic function with 1 method)

julia> @time f1 = fateman1(0);
  0.257372 seconds (210.58 k allocations: 9.605 MB, 2.08% gc time)

julia> @time f1 = fateman1(20);
  2.510551 seconds (3.43 k allocations: 20.112 MB, 0.18% gc time)

```

Another implementation of the same:

```julia
julia> function fateman2(degree::Int)
           T = Int128
           oneH = HomogeneousPolynomial(one(T), 0)
           s = TaylorN( [oneH, HomogeneousPolynomial([one(T),one(T),one(T),one(T)],1)], degree )
           s = s^degree
           # s is converted to order 2*ndeg
           s = TaylorN(s, 2*degree)
           return s^2 + s
       end
fateman2 (generic function with 1 method)

julia> @time f2 = fateman2(0);
  0.009598 seconds (7.26 k allocations: 320.872 KB)

julia> @time f2 = fateman2(20);
  1.219438 seconds (3.42 k allocations: 20.111 MB, 0.56% gc time)

julia> get_coeff(f2,[1,6,7,20]) # coef x^1 y^6 z^7 w^{20}
128358585324486316800

julia> sum(TaylorSeries.size_table) # number of distinct monomials
135751

```

The tests above show the necessity of using `Int128`, that
`fateman2` is nearly twice as fast as `fateman1`, and that the series has 135751
monomials on 4 variables.

Mathematica (version 10.2) requires 3.139831 seconds. Then,
with TaylorSeries v0.1.2, our implementation of `fateman1` is about 20% faster,
and `fateman2` is more than a factor 2 faster. (The original test by Fateman
corresponds to `fateman1`, which avoids optimizations in `^2`.)
