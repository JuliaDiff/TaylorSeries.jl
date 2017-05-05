# Background

---

## Introduction

[TaylorSeries.jl](https://github.com/lbenet/TaylorSeries.jl) is an implementation
of high-order
[automatic differentiation](http://en.wikipedia.org/wiki/Automatic_differentiation),
as presented in the book by W. Tucker [[1]](@ref refs). The general
idea is the following.

The Taylor series expansion of an analytical function
``f(t)`` with *one* independent variable ``t`` around ``t_0`` can be written as
```math
\begin{equation}
f(t) = f_0 + f_1 (t-t_0) + f_2 (t-t_0)^2 + \cdots + f_k (t-t_0)^k + \cdots,
\end{equation}
```
where ``f_0=f(t_0)``, and the Taylor coefficients ``f_k = f_k(t_0)`` are the
``k``-th *normalized derivatives* at ``t_0``:
```math
\begin{equation}
f_k = \frac{1}{k!} \frac{{\rm d}^k f} {{\rm d} t^k}(t_0).
\end{equation}
```
Thus, computing the high-order derivatives of ``f(t)`` is equivalent to computing
its Taylor expansion.

In the case of *many* independent variables the same statements hold, though
things become more subtle. Following Alex Haro's approach
[[2]](@ref refs), the Taylor
expansion is an infinite sum of *homogeneous polynomials* in the ``d`` independent
variables ``x_1, x_2, \dots, x_d``, which takes the form

```math
\begin{equation}
f_k (\mathbf{x_0}) = \sum_{m_1+\cdots+m_d = k}\, f_{m_1,\dots,m_d} \;\,
(x_1-x_{0_1})^{m_1} \cdots (x_d-x_{0_d})^{m_d} =
\sum_{|\mathbf{m}|=k} f_\mathbf{m}\, (\mathbf{x}-\mathbf{x_0})^\mathbf{m}.
\end{equation}
```

Here, ``\mathbf{m}\in \mathbb{N}^d`` is a multi-index of the ``k``-th order
homogeneous polynomial and ``\mathbf{x}=(x_1,x_2,\ldots,x_d)`` are the
``d`` independent variables.

In both cases, a Taylor series expansion can be represented by a
vector containing its coefficients. The difference between the cases of
one or more independent variables is that the
coefficients are real or complex numbers in the former case, but
homogeneous polynomials in the latter case. This motivates
the construction of the [`Taylor1`](@ref) and [`TaylorN`](@ref) types.

## Arithmetic operations

Arithmetic operations involving Taylor series can be expressed as
operations on the coefficients:

```math
\begin{eqnarray}
\label{eq:arith1}
(f(x) \pm g(x))_k & = & f_k \pm g_k , \\
\label{eq:arith2}
(f(x) \cdot g(x))_k & = & \sum_{i=0}^k f_i \, g_{k-i} , \\
\label{eq:arith3}
\Big( \frac{f(x)}{g(x)} \Big)_k & = & \frac{1}{g_0} \Big[ f_k -
\sum_{i=0}^{k-1} \big(\frac{f(x)}{g(x)}\big)_i \, g_{k-i} \Big]. \\
\end{eqnarray}
```

Equations (\ref{eq:arith1}-\ref{eq:arith3}) corresponds to a convolution.

## Elementary functions of polynomials

Consider a function ``y(t)`` that satisfies the ordinary differential equation
``\dot{y} = f(y)``, ``y(t_0)=y_0``, where ``t`` is the independent variable.
Writing ``y(t)`` and ``f(t)`` as Taylor polynomials of ``t``, substituting these in the
differential equation and equating equal powers of
the independent variable leads to the recursion relation

```math
\begin{equation}
\label{eq:rec}
y_{n+1} = \frac{f_n}{n+1}.
\end{equation}
```

Equation (\ref{eq:rec}) and the corresponding initial condition
``y(t_0)=y_0`` define a recurrence relation
for the Taylor coefficients of ``y(t)`` around ``t_0``.

The following are  examples of such recurrence relations for some
elementary functions:

```math
\begin{eqnarray}
p(t)=(f(t))^\alpha , \qquad &
  p_k & = \frac{1}{k \, f_0}\sum_{j=0}^{k-1}\big( \alpha(k-j)-j\big)
  \, f_{k-j} \, p_j; \\
e(t) = \exp(f(t)) , \qquad &
  e_k & = \frac{1}{k}\sum_{j=0}^{k-1} (k-j) \, f_{k-j} \, e_j; \\
l(t) = \log(f(t)) , \qquad &
  l_k & = \frac{1}{f_0}\big( f_k - \frac{1}{k}\sum_{j=1}^{k-1} j
    \, f_{k-j} \, l_j \big); \\
s(t) = \sin(f(t)) , \qquad &
  s_k & = \frac{1}{k}\sum_{j=0}^{k-1} (k-j) \, f_{k-j} \, c_j; \\
c(t) = \cos(f(t)) , \qquad &
  c_k & = -\frac{1}{k}\sum_{j=0}^{k-1} (k-j) \, f_{k-j} \, s_j.
\end{eqnarray}
```

The recursion relations for ``s(t) = \sin\big(f(t)\big)`` and
``c(t) = \cos\big(f(t)\big)`` depend
on each other; this reflects the fact that they are solutions of a second-order
differential equation.

All these relations hold for Taylor expansions in one
and more independent variables; in the latter case, the Taylor coefficients
``f_k`` are homogeneous polynomials of degree ``k``;
see [[2]](@ref refs).

## [References](@id refs)

[1] W. Tucker, *Validated Numerics: A Short Introduction to Rigorous
Computations*, Princeton University Press (2011).

[2] A. Haro, *Automatic differentiation methods in computational dynamical
systems: Invariant manifolds and normal forms of vector fields at fixed points*,
[preprint](http://www.maia.ub.es/~alex/admcds/admcds.pdf).
