# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#


# Functions
for T in (:Taylor1, :TaylorN)
    @eval begin
        function exp(a::$T)
            order = a.order
            aux = exp(constant_term(a))
            aa = one(aux) * a
            c = $T( aux, order )
            for k in eachindex(a)
                exp!(c, aa, k)
            end
            return c
        end

        function log(a::$T)
            iszero(constant_term(a)) && throw(DomainError(a, 
                    """The 0-th order coefficient must be non-zero in order to expand `log` around 0."""))

            order = a.order
            aux = log(constant_term(a))
            aa = one(aux) * a
            c = $T( aux, order )
            for k in eachindex(a)
                log!(c, aa, k)
            end
            return c
        end

        sin(a::$T) = sincos(a)[1]
        cos(a::$T) = sincos(a)[2]
        function sincos(a::$T)
            order = a.order
            aux = sin(constant_term(a))
            aa = one(aux) * a
            s = $T( aux, order )
            c = $T( cos(constant_term(a)), order )
            for k in eachindex(a)
                sincos!(s, c, aa, k)
            end
            return s, c
        end

        function tan(a::$T)
            order = a.order
            aux = tan(constant_term(a))
            aa = one(aux) * a
            c = $T(aux, order)
            c2 = $T(aux^2, order)
            for k in eachindex(a)
                tan!(c, aa, c2, k)
            end
            return c
        end

        function asin(a::$T)
            a0 = constant_term(a)
            a0^2 == one(a0) && throw(DomainError(a, 
                """Series expansion of asin(x) diverges at x = ±1."""))

            order = a.order
            aux = asin(a0)
            aa = one(aux) * a
            c = $T( aux, order )
            r = $T( sqrt(1 - a0^2), order )
            for k in eachindex(a)
                asin!(c, aa, r, k)
            end
            return c
        end

        function acos(a::$T)
            a0 = constant_term(a)
            a0^2 == one(a0) && throw(DomainError(a, 
            """Series expansion of asin(x) diverges at x = ±1."""))

            order = a.order
            aux = acos(a0)
            aa = one(aux) * a
            c = $T( aux, order )
            r = $T( sqrt(1 - a0^2), order )
            for k in eachindex(a)
                acos!(c, aa, r, k)
            end
            return c
        end

        function atan(a::$T)
            order = a.order
            a0 = constant_term(a)
            aux = atan(a0)
            aa = one(aux) * a
            c = $T( aux, order)
            r = $T(one(aux) + a0^2, order)
            iszero(constant_term(r)) && throw(DomainError(a, 
                """Series expansion of atan(x) diverges at x = ±im."""))

            for k in eachindex(a)
                atan!(c, aa, r, k)
            end
            return c
        end

        function atan(a::$T, b::$T)
            c = atan(a/b)
            c[0] = atan(constant_term(a), constant_term(b))
            return c
        end

        sinh(a::$T) = sinhcosh(a)[1]
        cosh(a::$T) = sinhcosh(a)[2]
        function sinhcosh(a::$T)
            order = a.order
            aux = sinh(constant_term(a))
            aa = one(aux) * a
            s = $T( aux, order)
            c = $T( cosh(constant_term(a)), order)
            for k in eachindex(a)
                sinhcosh!(s, c, aa, k)
            end
            return s, c
        end

        function tanh(a::$T)
            order = a.order
            aux = tanh( constant_term(a) )
            aa = one(aux) * a
            c = $T( aux, order)
            c2 = $T( aux^2, order)
            for k in eachindex(a)
                tanh!(c, aa, c2, k)
            end
            return c
        end


        function asinh(a::$T)
            order = a.order
            a0 = constant_term(a)
            aux = asinh(a0)
            aa = one(aux) * a
            c = $T( aux, order )
            r = $T( sqrt(a0^2 + 1), order )
            iszero(constant_term(r)) && throw(DomainError(a,
                """Series expansion of asinh(x) diverges at x = ±im."""))
            for k in eachindex(a)
                asinh!(c, aa, r, k)
            end
            return c
        end

        function acosh(a::$T)
            a0 = constant_term(a)
            a0^2 == one(a0) && throw(DomainError(a,
                """Series expansion of acosh(x) diverges at x = ±1."""))

            order = a.order
            aux = acosh(a0)
            aa = one(aux) * a
            c = $T( aux, order )
            r = $T( sqrt(a0^2 - 1), order )
            for k in eachindex(a)
                acosh!(c, aa, r, k)
            end
            return c
        end

        function atanh(a::$T)
            order = a.order
            a0 = constant_term(a)
            aux = atanh(a0)
            aa = one(aux) * a
            c = $T( aux, order)
            r = $T(one(aux) - a0^2, order)
            iszero(constant_term(r)) && throw(DomainError(a,
                """Series expansion of atanh(x) diverges at x = ±1."""))

            for k in eachindex(a)
                atanh!(c, aa, r, k)
            end
            return c
        end



    end
end


# Recursive functions (homogeneous coefficients)
for T in (:Taylor1, :TaylorN)
    @eval begin
        @inline function identity!(c::$T{T}, a::$T{T}, k::Int) where {T<:Number}
            @inbounds c[k] = identity(a[k])
            return nothing
        end

        @inline function zero!(c::$T{T}, a::$T{T}, k::Int) where {T<:Number}
            @inbounds c[k] = zero(a[k])
            return nothing
        end

        @inline function one!(c::$T{T}, a::$T{T}, k::Int) where {T<:Number}
            if k == 0
                @inbounds c[0] = one(a[0])
            else
                @inbounds c[k] = zero(a[k])
            end
            return nothing
        end

        @inline function abs!(c::$T{T}, a::$T{T}, k::Int) where {T<:Number}
            z = zero(constant_term(a))
            if constant_term(constant_term(a)) > constant_term(z)
                return add!(c, a, k)
            elseif constant_term(constant_term(a)) < constant_term(z)
                return subst!(c, a, k)
            else
                throw(DomainError(a,
                """The 0th order coefficient must be non-zero
                (abs(x) is not differentiable at x=0)."""))
            end
            return nothing
        end

        @inline abs2!(c::$T{T}, a::$T{T}, k::Int) where {T<:Number} = sqr!(c, a, k)

        @inline function exp!(c::$T{T}, a::$T{T}, k::Int) where {T<:Number}
            if k == 0
                @inbounds c[0] = exp(constant_term(a))
                return nothing
            end
            if $T == Taylor1
                @inbounds c[k] = k * a[k] * c[0]
            else
                @inbounds mul!(c[k], k * a[k], c[0])
            end
            @inbounds for i = 1:k-1
                if $T == Taylor1
                    c[k] += (k-i) * a[k-i] * c[i]
                else
                    mul!(c[k], (k-i) * a[k-i], c[i])
                end
            end
            @inbounds c[k] = c[k] / k
            return nothing
        end

        @inline function log!(c::$T{T}, a::$T{T}, k::Int) where {T<:Number}
            if k == 0
                @inbounds c[0] = log(constant_term(a))
                return nothing
            elseif k == 1
                @inbounds c[1] = a[1] / constant_term(a)
                return nothing
            end

            if $T == Taylor1
                @inbounds c[k] = (k-1) * a[1] * c[k-1]
            else
                @inbounds mul!(c[k], (k-1)*a[1], c[k-1])
            end
            @inbounds for i = 2:k-1
                if $T == Taylor1
                    c[k] += (k-i) * a[i] * c[k-i]
                else
                    mul!(c[k], (k-i)*a[i], c[k-i])
                end
            end
            @inbounds c[k] = (a[k] - c[k]/k) / constant_term(a)
            return nothing
        end

        @inline function sincos!(s::$T{T}, c::$T{T}, a::$T{T}, k::Int) where {T<:Number}
            if k == 0
                a0 = constant_term(a)
                @inbounds s[0], c[0] = sin( a0 ), cos( a0 )
                return nothing
            end

            x = a[1]
            if $T == Taylor1
                @inbounds s[k] = x * c[k-1]
                @inbounds c[k] = -x * s[k-1]
            else
                mul!(s[k], x, c[k-1])
                mul!(c[k], -x, s[k-1])
            end
            @inbounds for i = 2:k
                x = i * a[i]
                if $T == Taylor1
                    s[k] += x * c[k-i]
                    c[k] -= x * s[k-i]
                else
                    mul!(s[k], x, c[k-i])
                    mul!(c[k], -x, s[k-i])
                end
            end

            @inbounds s[k] = s[k] / k
            @inbounds c[k] = c[k] / k
            return nothing
        end

        @inline function tan!(c::$T{T}, a::$T{T}, c2::$T{T}, k::Int) where {T<:Number}
            if k == 0
                @inbounds aux = tan( constant_term(a) )
                @inbounds c[0] = aux
                @inbounds c2[0] = aux^2
                return nothing
            end

            if $T == Taylor1
                @inbounds c[k] = (k-1) * a[k-1] * c2[1]
            else
                @inbounds mul!(c[k], (k-1) * a[k-1], c2[1])
            end
            @inbounds for i = 1:k-1
                if $T == Taylor1
                    c[k] += (k-i) * a[k-i] * c2[i]
                else
                    mul!(c[k], (k-i) * a[k-i], c2[i])
                end
            end
            @inbounds c[k] = a[k] + c[k]/k
            sqr!(c2, c, k)

            return nothing
        end

        @inline function asin!(c::$T{T}, a::$T{T}, r::$T{T}, k::Int) where {T<:Number}
            if k == 0
                a0 = constant_term(a)
                @inbounds c[0] = asin( a0 )
                @inbounds r[0] = sqrt( 1 - a0^2 )
                return nothing
            end

            if $T == Taylor1
                @inbounds c[k] = (k-1) * r[1] * c[k-1]
            else
                @inbounds mul!(c[k], (k-1) * r[1], c[k-1])
            end
            @inbounds for i in 2:k-1
                if $T == Taylor1
                    c[k] += (k-i) * r[i] * c[k-i]
                else
                    mul!(c[k], (k-i) * r[i], c[k-i])
                end
            end
            sqrt!(r, 1-a^2, k)
            @inbounds c[k] = (a[k] - c[k]/k) / constant_term(r)
            return nothing
        end

        @inline function acos!(c::$T{T}, a::$T{T}, r::$T{T}, k::Int) where {T<:Number}
            if k == 0
                a0 = constant_term(a)
                @inbounds c[0] = acos( a0 )
                @inbounds r[0] = sqrt( 1 - a0^2 )
                return nothing
            end

            if $T == Taylor1
                @inbounds c[k] = (k-1) * r[1] * c[k-1]
            else
                @inbounds mul!(c[k], (k-1) * r[1], c[k-1])
            end
            @inbounds for i in 2:k-1
                if $T == Taylor1
                    c[k] += (k-i) * r[i] * c[k-i]
                else
                    mul!(c[k], (k-i) * r[i], c[k-i])
                end
            end
            sqrt!(r, 1-a^2, k)
            @inbounds c[k] = -(a[k] + c[k]/k) / constant_term(r)
            return nothing
        end

        @inline function atan!(c::$T{T}, a::$T{T}, r::$T{T}, k::Int) where {T<:Number}
            if k == 0
                a0 = constant_term(a)
                @inbounds c[0] = atan( a0 )
                @inbounds r[0] = 1 + a0^2
                return nothing
            end

            if $T == Taylor1
                @inbounds c[k] = (k-1) * r[1] * c[k-1]
            else
                @inbounds mul!(c[k], (k-1) * r[1], c[k-1])
            end
            @inbounds for i in 2:k-1
                if $T == Taylor1
                    c[k] += (k-i) * r[i] * c[k-i]
                else
                    mul!(c[k], (k-i) * r[i], c[k-i])
                end
            end
            @inbounds sqr!(r, a, k)
            @inbounds c[k] = (a[k] - c[k]/k) / constant_term(r)
            return nothing
        end

        @inline function sinhcosh!(s::$T{T}, c::$T{T}, a::$T{T}, k::Int) where {T<:Number}
            if k == 0
                @inbounds s[0] = sinh( constant_term(a) )
                @inbounds c[0] = cosh( constant_term(a) )
                return nothing
            end

            x = a[1]
            if $T == Taylor1
                @inbounds s[k] = x * c[k-1]
                @inbounds c[k] = x * s[k-1]
            else
                @inbounds mul!(s[k], x, c[k-1])
                @inbounds mul!(c[k], x, s[k-1])
            end
            @inbounds for i = 2:k
                x = i * a[i]
                if $T == Taylor1
                    s[k] += x * c[k-i]
                    c[k] += x * s[k-i]
                else
                    mul!(s[k], x, c[k-i])
                    mul!(c[k], x, s[k-i])
                end
            end
            s[k] = s[k] / k
            c[k] = c[k] / k
            return nothing
        end

        @inline function tanh!(c::$T{T}, a::$T{T}, c2::$T{T}, k::Int) where {T<:Number}
            if k == 0
                @inbounds aux = tanh( constant_term(a) )
                @inbounds c[0] = aux
                @inbounds c2[0] = aux^2
                return nothing
            end

            if $T == Taylor1
                @inbounds c[k] = k * a[k] * c2[0]
            else
                @inbounds mul!(c[k], k * a[k], c2[0])
            end
            @inbounds for i = 1:k-1
                if $T == Taylor1
                    c[k] += (k-i) * a[k-i] * c2[i]
                else
                    mul!(c[k], (k-i) * a[k-i], c2[i])
                end
            end
            @inbounds c[k] = a[k] - c[k]/k
            sqr!(c2, c, k)

            return nothing
        end

        @inline function asinh!(c::$T{T}, a::$T{T}, r::$T{T}, k::Int) where {T<:Number}
            if k == 0
                a0 = constant_term(a)
                @inbounds c[0] = asinh( a0 )
                @inbounds r[0] = sqrt( a0^2 + 1 )
                return nothing
            end

            if $T == Taylor1
                @inbounds c[k] = (k-1) * r[1] * c[k-1]
            else
                @inbounds mul!(c[k], (k-1) * r[1], c[k-1])
            end
            @inbounds for i in 2:k-1
                if $T == Taylor1
                    c[k] += (k-i) * r[i] * c[k-i]
                else
                    mul!(c[k], (k-i) * r[i], c[k-i])
                end
            end
            sqrt!(r, a^2+1, k)
            @inbounds c[k] = (a[k] - c[k]/k) / constant_term(r)
            return nothing
        end

        @inline function acosh!(c::$T{T}, a::$T{T}, r::$T{T}, k::Int) where {T<:Number}
            if k == 0
                a0 = constant_term(a)
                @inbounds c[0] = acosh( a0 )
                @inbounds r[0] = sqrt( a0^2 - 1 )
                return nothing
            end

            if $T == Taylor1
                @inbounds c[k] = (k-1) * r[1] * c[k-1]
            else
                @inbounds mul!(c[k], (k-1) * r[1], c[k-1])
            end
            @inbounds for i in 2:k-1
                if $T == Taylor1
                    c[k] += (k-i) * r[i] * c[k-i]
                else
                    mul!(c[k], (k-i) * r[i], c[k-i])
                end
            end
            sqrt!(r, a^2-1, k)
            @inbounds c[k] = (a[k] - c[k]/k) / constant_term(r)
            return nothing
        end

        @inline function atanh!(c::$T{T}, a::$T{T}, r::$T{T}, k::Int) where {T<:Number}
            if k == 0
                a0 = constant_term(a)
                @inbounds c[0] = atanh( a0 )
                @inbounds r[0] = 1 - a0^2
                return nothing
            end

            if $T == Taylor1
                @inbounds c[k] = (k-1) * r[1] * c[k-1]
            else
                @inbounds mul!(c[k], (k-1) * r[1], c[k-1])
            end
            @inbounds for i in 2:k-1
                if $T == Taylor1
                    c[k] += (k-i) * r[i] * c[k-i]
                else
                    mul!(c[k], (k-i) * r[i], c[k-i])
                end
            end
            @inbounds sqr!(r, a, k)
            @inbounds c[k] = (a[k] + c[k]/k) / constant_term(r)
            return nothing
        end


    end
end


@doc doc"""
    inverse(f)

Return the Taylor expansion of ``f^{-1}(t)``, of order `N = f.order`,
for `f::Taylor1` polynomial if the first coefficient of `f` is zero.
Otherwise, a `DomainError` is thrown.

The algorithm implements Lagrange inversion at ``t=0`` if ``f(0)=0``:
```math
\begin{equation*}
f^{-1}(t) = \sum_{n=1}^{N} \frac{t^n}{n!} \left.
    \frac{{\rm d}^{n-1}}{{\rm d} z^{n-1}}\left(\frac{z}{f(z)}\right)^n
    \right\vert_{z=0}.
\end{equation*}
```

""" inverse

function inverse(f::Taylor1{T}) where {T<:Number}
    if f[0] != zero(T)
        throw(DomainError(f,
        """
        Evaluation of Taylor1 series at 0 is non-zero. For high accuracy, revert
        a Taylor1 series with first coefficient 0 and re-expand about f(0).
        """))
    end
    z = Taylor1(T,f.order)
    zdivf = z/f
    zdivfpown = zdivf
    S = eltype(zdivf)
    coeffs = zeros(S,f.order+1)

    @inbounds for n in 1:f.order
        coeffs[n+1] = zdivfpown[n-1]/n
        zdivfpown *= zdivf
    end
    Taylor1(coeffs, f.order)
end



# Documentation for the recursion relations
@doc doc"""
    exp!(c, a, k) --> nothing

Update the `k-th` expansion coefficient `c[k+1]` of `c = exp(a)`
for both `c` and `a` either `Taylor1` or `TaylorN`.

The coefficients are given by

```math
\begin{equation*}
c_k = \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} c_j.
\end{equation*}
```

""" exp!


@doc doc"""
    log!(c, a, k) --> nothing

Update the `k-th` expansion coefficient `c[k+1]` of `c = log(a)`
for both `c` and `a` either `Taylor1` or `TaylorN`.

The coefficients are given by

```math
\begin{equation*}
c_k = \frac{1}{a_0} \big(a_k - \frac{1}{k} \sum_{j=0}^{k-1} j a_{k-j} c_j \big).
\end{equation*}
```

""" log!


@doc doc"""
    sincos!(s, c, a, k) --> nothing

Update the `k-th` expansion coefficients `s[k+1]` and `c[k+1]`
of `s = sin(a)` and `c = cos(a)` simultaneously, for `s`, `c` and `a`
either `Taylor1` or `TaylorN`.

The coefficients are given by

```math
\begin{eqnarray*}
s_k &=&  \frac{1}{k}\sum_{j=0}^{k-1} (k-j) a_{k-j} c_j ,\\
c_k &=& -\frac{1}{k}\sum_{j=0}^{k-1} (k-j) a_{k-j} s_j.
\end{eqnarray*}
```

""" sincos!


@doc doc"""
    tan!(c, a, p, k::Int) --> nothing

Update the `k-th` expansion coefficients `c[k+1]` of `c = tan(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `p = c^2` and
is passed as an argument for efficiency.

The coefficients are given by

```math
\begin{equation*}
c_k = a_k + \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} p_j.
\end{equation*}
```

""" tan!


@doc doc"""
    asin!(c, a, r, k)

Update the `k-th` expansion coefficients `c[k+1]` of `c = asin(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `r = sqrt(1-c^2)` and
is passed as an argument for efficiency.

```math
\begin{equation*}
c_k = \frac{1}{ \sqrt{r_0} }
    \big( a_k - \frac{1}{k} \sum_{j=1}^{k-1} j r_{k-j} c_j \big).
\end{equation*}
```

""" asin!


@doc doc"""
    acos!(c, a, r, k)

Update the `k-th` expansion coefficients `c[k+1]` of `c = acos(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `r = sqrt(1-c^2)` and
is passed as an argument for efficiency.

```math
\begin{equation*}
c_k = - \frac{1}{ r_0 }
    \big( a_k - \frac{1}{k} \sum_{j=1}^{k-1} j r_{k-j} c_j \big).
\end{equation*}
```

""" acos!


@doc doc"""
    atan!(c, a, r, k)

Update the `k-th` expansion coefficients `c[k+1]` of `c = atan(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `r = 1+a^2` and
is passed as an argument for efficiency.

```math
\begin{equation*}
c_k = \frac{1}{r_0}\big(a_k - \frac{1}{k} \sum_{j=1}^{k-1} j r_{k-j} c_j\big).
\end{equation*}
```

""" atan!


@doc doc"""
    sinhcosh!(s, c, a, k)

Update the `k-th` expansion coefficients `s[k+1]` and `c[k+1]`
of `s = sinh(a)` and `c = cosh(a)` simultaneously, for `s`, `c` and `a`
either `Taylor1` or `TaylorN`.

The coefficients are given by

```math
\begin{eqnarray*}
s_k &=& \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} c_j, \\
c_k &=& \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} s_j.
\end{eqnarray*}
```

""" sinhcosh!


@doc doc"""
    tanh!(c, a, p, k)

Update the `k-th` expansion coefficients `c[k+1]` of `c = tanh(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `p = a^2` and
is passed as an argument for efficiency.

```math
\begin{equation*}
c_k = a_k - \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} p_j.
\end{equation*}
```

""" tanh!



@doc doc"""
    asinh!(c, a, r, k)

Update the `k-th` expansion coefficients `c[k+1]` of `c = asinh(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `r = sqrt(1-c^2)` and
is passed as an argument for efficiency.

```math
\begin{equation*}
c_k = \frac{1}{ \sqrt{r_0} }
    \big( a_k - \frac{1}{k} \sum_{j=1}^{k-1} j r_{k-j} c_j \big).
\end{equation*}
```

""" asinh!


@doc doc"""
    acosh!(c, a, r, k)

Update the `k-th` expansion coefficients `c[k+1]` of `c = acosh(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `r = sqrt(c^2-1)` and
is passed as an argument for efficiency.

```math
\begin{equation*}
c_k = \frac{1}{ r_0 }
    \big( a_k - \frac{1}{k} \sum_{j=1}^{k-1} j r_{k-j} c_j \big).
\end{equation*}
```

""" acosh!


@doc doc"""
    atanh!(c, a, r, k)

Update the `k-th` expansion coefficients `c[k+1]` of `c = atanh(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `r = 1-a^2` and
is passed as an argument for efficiency.

```math
\begin{equation*}
c_k = \frac{1}{r_0}\big(a_k + \frac{1}{k} \sum_{j=1}^{k-1} j r_{k-j} c_j\big).
\end{equation*}
```

""" atanh!
