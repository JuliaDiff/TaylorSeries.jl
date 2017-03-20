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
        ## exp ##
        function exp(a::$T)
            order = max_order(a)
            @inbounds aux = exp( constant_term(a) )
            c = $T( aux, order )
            @inbounds for k = 1:order
                exp!(c, a, k)
            end
            return c
        end

        ## log ##
        function log(a::$T)
            constant_term(a) == zero(constant_term(a)) &&
                throw(ArgumentError("""
                    The 0-th order `TaylorN` coefficient must be non-zero
                    in order to expand `log` around 0.
                    """))

            order = max_order(a)
            @inbounds aux = log( constant_term(a) )
            c = $T( aux, order )
            @inbounds for k = 1:order
                log!(c, a, k)
            end
            return c
        end

        ## sin and cos ##
        sin(a::$T) = sincos(a)[1]
        cos(a::$T) = sincos(a)[2]
        function sincos(a::$T)
            order = max_order(a)
            @inbounds s = $T( sin(constant_term(a)), order )
            @inbounds c = $T( cos(constant_term(a)), order )
            @inbounds for k = 1:order
                sincos!(s, c, a, k)
            end
            return s, c
        end

        ## tan ##
        function tan(a::$T)
            order = max_order(a)
            @inbounds aux = tan(constant_term(a))
            c = $T(aux, order)
            c2 = $T(aux^2, order)
            @inbounds for k = 1:order
                tan!(c, a, c2, k)
            end
            return c
        end

        ## asin ##
        function asin(a::$T)
            a0 = constant_term(a)
            a0^2 == one(a0) && throw(ArgumentError(
                """
                Recursion formula diverges due to vanishing `sqrt`
                in the denominator.
                """))

            order = max_order(a)
            c = $T( asin(a0), order )
            r = $T( sqrt(1 - a0^2), order )
            @inbounds for k in 1:order
                asin!(c, a, r, k)
            end
            return c
        end

        ## acos ##
        function acos(a::$T)
            a0 = constant_term(a)
            a0^2 == one(a0) && throw(ArgumentError(
                """
                Recursion formula diverges due to vanishing `sqrt`
                in the denominator.
                """))

            order = max_order(a)
            c = $T( asin(a0), order )
            r = $T( sqrt(1 - a0^2), order )
            @inbounds for k in 1:order
                acos!(c, a, r, k)
            end
            @inbounds c[1] = acos(a0)
            return c
        end

        ## atan ##
        function atan(a::$T)
            order = max_order(a)
            a0 = constant_term(a)
            c = $T( atan(a0), order)
            r = $T(1 + a0^2, order)
            constant_term(r) == zero(constant_term(a)) &&
                throw(ArgumentError(
                    """
                    Recursion formula has a pole.
                    """))

            @inbounds for k in 1:order
                atan!(c, a, r, k)
            end
            return c
        end

        ## sinh and cosh ##
        sinh(a::$T) = sinhcosh(a)[1]
        cosh(a::$T) = sinhcosh(a)[2]
        function sinhcosh(a::$T)
            order = max_order(a)
            a0 = constant_term(a)
            s = $T( sinh(a0), order)
            c = $T( cosh(a0), order)
            @inbounds for k = 1:order
                sinhcosh!(s, c, a, k)
            end
            return s, c
        end

        ## tanh ##
        function tanh(a::$T)
            order = max_order(a)
            @inbounds aux = tanh( constant_term(a) )
            c = $T( aux, order)
            c2 = $T( aux^2, order)
            @inbounds for k = 1:order
                tanh!(c, a, c2, k)
            end
            return c
        end
    end
end


# Recursive functions (homogeneous coefficients)
for T in (:Taylor1, :TaylorN)
    @eval begin
        function exp!(c::$T, a::$T, k::Int)
            if $T == Taylor1
                if k == 0
                    @inbounds c[1] = exp(a[1])
                    return nothing
                end
                @inbounds c[k+1] = zero(c[1])
            else
                if k == 0
                    @inbounds c[1] = exp(a[1][1])
                    return nothing
                end
                @inbounds c[k+1] = HomogeneousPolynomial(zero(c[1][1]), k)
            end

            @inbounds for i = 0:k-1
                c[k+1] += (k-i) * a[k-i+1] * c[i+1]
            end
            @inbounds c[k+1] = c[k+1] / k

            return nothing
        end

        function log!(c::$T, a::$T, k::Int)
            if $T == Taylor1
                if k == 0
                    @inbounds c[1] = log(a[1])
                    return nothing
                end
                @inbounds c[k+1] = zero(c[1])
            else
                if k == 0
                    @inbounds c[1] = log(a[1][1])
                    return nothing
                end
                @inbounds c[k+1] = HomogeneousPolynomial(zero(c[1][1]), k)
            end

            @inbounds for i = 1:k-1
                c[k+1] += (k-i) * a[i+1] * c[k-i+1]
            end

            if $T == Taylor1
                @inbounds c[k+1] = (a[k+1] -c[k+1]/k) / a[1]
            else
                @inbounds c[k+1] = (a[k+1] -c[k+1]/k) / a[1][1]
            end
            return nothing
        end

        function sincos!(s::$T, c::$T, a::$T, k::Int)
            if $T == Taylor1
                if k == 0
                    @inbounds s[1], c[1] = sin( a[1] ), cos( a[1] )
                    return nothing
                end
                @inbounds s[k+1] = zero(s[1])
                @inbounds c[k+1] = zero(c[1])
            else
                if k == 0
                    @inbounds s[1], c[1] = sin( a[1][1] ), cos( a[1][1] )
                    return nothing
                end
                @inbounds s[k+1] = HomogeneousPolynomial(zero(s[1][1]), k)
                @inbounds c[k+1] = HomogeneousPolynomial(zero(c[1][1]), k)
            end

            @inbounds for i = 1:k
                x = i * a[i+1]
                s[k+1] += x * c[k-i+1]
                c[k+1] -= x * s[k-i+1]
            end

            @inbounds s[k+1] = s[k+1]/k
            @inbounds c[k+1] = c[k+1]/k
            return nothing
        end

        function tan!(c::$T, a::$T, c2::$T, k::Int)
            if $T == Taylor1
                if k == 0
                    @inbounds aux = tan( a[1] )
                    @inbounds c[1] = aux
                    @inbounds c2[1] = aux^2
                    return nothing
                end
                @inbounds c[k+1] = zero(c[1])
            else
                if k == 0
                    @inbounds aux = tan( a[1][1] )
                    @inbounds c[1] = aux
                    @inbounds c2[1] = aux^2
                    return nothing
                end
                @inbounds c[k+1] = HomogeneousPolynomial(zero(c[1][1]), k)
            end

            @inbounds for i = 0:k-1
                c[k+1] += (k-i)*a[k-i+1]*c2[i+1]
            end
            @inbounds c[k+1] = a[k+1] + c[k+1]/k
            sqr!(c2, c, k)

            return nothing
        end

        function asin!(c::$T, a::$T, r::$T, k::Int)
            if $T == Taylor1
                if k == 0
                    @inbounds c[1] = asin( a[1] )
                    @inbounds r[1] = sqrt( 1 - a[1]^2)
                    return nothing
                end
                @inbounds c[k+1] = zero(c[1])
            else
                if k == 0
                    @inbounds c[1] = asin( a[1][1] )
                    @inbounds r[1] = sqrt( 1 - a[1][1]^2)
                    return nothing
                end
                @inbounds c[k+1] = HomogeneousPolynomial(zero(c[1][1]), k)
            end

            @inbounds for i in 1:k-1
                c[k+1] += (k-i) * r[i+1] * c[k-i+1]
            end
            sqrt!(r, 1-a^2, k)
            if $T == Taylor1
                @inbounds c[k+1] = (a[k+1] - c[k+1]/k) / r[1]
            else
                @inbounds c[k+1] = (a[k+1] - c[k+1]/k) / r[1][1]
            end
            return nothing
        end

        function acos!(c::$T, a::$T, r::$T, k::Int)
            if $T == Taylor1
                if k == 0
                    @inbounds c[1] = acos( a[1] )
                    @inbounds r[1] = sqrt( 1 - a[1]^2)
                    return nothing
                end
                @inbounds c[k+1] = zero(c[1])
            else
                if k == 0
                    @inbounds c[1] = acos( a[1][1] )
                    @inbounds r[1] = sqrt( 1 - a[1][1]^2)
                    return nothing
                end
                @inbounds c[k+1] = HomogeneousPolynomial(zero(c[1][1]), k)
            end

            asin!(c, -a, r, k)
            return nothing
        end

        function atan!(c::$T, a::$T, r::$T, k::Int)
            if $T == Taylor1
                if k == 0
                    @inbounds c[1] = atan( a[1] )
                    @inbounds r[1] = 1 + a[1]^2
                    return nothing
                end
                @inbounds c[k+1] = zero(c[1])
            else
                if k == 0
                    @inbounds c[1] = atan( a[1][1] )
                    @inbounds r[1] = 1 + a[1][1]^2
                    return nothing
                end
                @inbounds c[k+1] = HomogeneousPolynomial(zero(c[1][1]), k)
            end

            @inbounds for i in 1:k-1
                c[k+1] += (k-i) * r[i+1] * c[k-i+1]
            end
            @inbounds sqr!(r, a, k)
            if $T == Taylor1
                @inbounds c[k+1] = (a[k+1] - c[k+1]/k) / r[1]
            else
                @inbounds c[k+1] = (a[k+1] - c[k+1]/k) / r[1][1]
            end

            return nothing
        end

        function sinhcosh!(s::$T, c::$T, a::$T, k::Int)
            if $T == Taylor1
                if k == 0
                    @inbounds s[1] = asinh( a[1] )
                    @inbounds c[1] = acosh( a[1] )
                    return nothing
                end
                @inbounds s[k+1] = zero(s[1])
                @inbounds c[k+1] = zero(c[1])
            else
                if k == 0
                    @inbounds s[1] = asinh( a[1][1] )
                    @inbounds c[1] = acosh( a[1][1] )
                    return nothing
                end
                @inbounds s[k+1] = HomogeneousPolynomial(zero(s[1][1]), k)
                @inbounds c[k+1] = HomogeneousPolynomial(zero(c[1][1]), k)
            end

            @inbounds for i = 1:k
                x = i * a[i+1]
                s[k+1] += x * c[k-i+1]
                c[k+1] += x * s[k-i+1]
            end
            s[k+1] = s[k+1]/k
            c[k+1] = c[k+1]/k
            return nothing
        end

        function tanh!(c::$T, a::$T, c2::$T, k::Int)
            if $T == Taylor1
                if k == 0
                    @inbounds aux = tanh( a[1] )
                    @inbounds c[1] = aux
                    @inbounds c2[1] = aux^2
                    return nothing
                end
                @inbounds c[k+1] = zero(c[1])
            else
                if k == 0
                    @inbounds aux = tanh( a[1][1] )
                    @inbounds c[1] = aux
                    @inbounds c2[1] = aux^2
                    return nothing
                end
                @inbounds c[k+1] = HomogeneousPolynomial(zero(c[1][1]), k)
            end

            @inbounds for i = 0:k-1
                c[k+1] += (k-i)*a[k-i+1]*c2[i+1]
            end
            @inbounds c[k+1] = a[k+1] - c[k+1]/k
            sqr!(c2, c, k)

            return nothing
        end
    end
end


doc"""
    reverse(f)

Return the Taylor expansion of $f^{-1}(t)$, of order `N = f.order`,
for `f::Taylor1` polynomial if the first coefficient of `f` is zero.
Otherwise, an `ArgumentError` is thrown.

The algorithm implements Lagrange inversion at $t=0$ if $f(0)=0$:
\begin{equation*}
f^{-1}(t) = \sum_{n=1}^{N} \frac{t^n}{n!} \left.
    \frac{{\rm d}^{n-1}}{{\rm d} z^{n-1}}\left(\frac{z}{f(z)}\right)^n
    \right\vert_{z=0}.
\end{equation*}
"""
function reverse{T<:Number}(f::Taylor1{T})
    if f[1] != zero(T)
        throw(ArgumentError(
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

    coeffs[1] = zero(S)
    @inbounds for n = 1:f.order
        coeffs[n+1] = zdivfpown[n]/n
        zdivfpown *= zdivf
    end
    Taylor1(coeffs, f.order)
end



# Documentation for the recursion relations
doc"""
    exp!(c, a, k) --> nothing

Update the `k-th` expansion coefficient `c[k+1]` of `c = exp(a)`
for both `c` and `a` either `Taylor1` or `TaylorN`.

The coefficients are given by

\begin{equation*}
c_k = \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} c_j.
\end{equation*}

""" exp!


doc"""
    log!(c, a, k) --> nothing

Update the `k-th` expansion coefficient `c[k+1]` of `c = log(a)`
for both `c` and `a` either `Taylor1` or `TaylorN`.

The coefficients are given by

\begin{equation*}
c_k = \frac{1}{a_0} (a_k - \frac{1}{k} \sum_{j=0}^{k-1} j a_{k-j} c_j ).
\end{equation*}

""" log!


doc"""
    sincos!(s, c, a, k) --> nothing

Update the `k-th` expansion coefficients `s[k+1]` and `c[k+1]`
of `s = sin(a)` and `c = cos(a)` simultaneously, for `s`, `c` and `a`
either `Taylor1` or `TaylorN`.

The coefficients are given by

\\begin{eqnarray*}
s_k &=&  \\frac{1}{k}\\sum_{j=0}^{k-1} (k-j) a_{k-j} c_j ,\\\\
c_k &=& -\\frac{1}{k}\\sum_{j=0}^{k-1} (k-j) a_{k-j} s_j.
\\end{eqnarray*}

""" sincos!


doc"""
    tan!(c, a, p, k::Int) --> nothing

Update the `k-th` expansion coefficients `c[k+1]` of `c = tan(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `p = c^2` and
is passed as an argument for efficiency.

The coefficients are given by

\begin{equation*}
c_k = a_k + \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} p_j.
\end{equation*}

""" tan!


doc"""
    asin!(c, a, r, k)

Update the `k-th` expansion coefficients `c[k+1]` of `c = asin(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `r = sqrt(1-c^2)` and
is passed as an argument for efficiency.

\begin{equation*}
c_k = \frac{1}{ \sqrt{r_0} }
    \big( a_k - \frac{1}{k} \sum_{j=1}^{k-1} j r_{k-j} c_j \big).
\end{equation*}

""" asin!


doc"""
    acos!(c, a, r, k)

Update the `k-th` expansion coefficients `c[k+1]` of `c = acos(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `r = sqrt(1-c^2)` and
is passed as an argument for efficiency.


\begin{equation*}
c_k = - \frac{1}{ r_0 }
    \big( a_k - \frac{1}{k} \sum_{j=1}^{k-1} j r_{k-j} c_j \big).
\end{equation*}

""" acos!


doc"""
    atan!(c, a, r, k)

Update the `k-th` expansion coefficients `c[k+1]` of `c = atan(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `r = 1+a^2` and
is passed as an argument for efficiency.

\begin{equation*}
c_k = \frac{1}{r_0}(a_k - \frac{1}{k} \sum_{j=1}^{k-1} j r_{k-j} c_j).
\end{equation*}

""" atan!


doc"""
    sinhcosh!(s, c, a, k)

Update the `k-th` expansion coefficients `s[k+1]` and `c[k+1]`
of `s = sinh(a)` and `c = cosh(a)` simultaneously, for `s`, `c` and `a`
either `Taylor1` or `TaylorN`.

The coefficients are given by

\\begin{eqnarray*}
s_k &=& \\frac{1}{k} \\sum_{j=0}^{k-1} (k-j) a_{k-j} c_j, \\\\
c_k &=& \\frac{1}{k} \\sum_{j=0}^{k-1} (k-j) a_{k-j} s_j.
\\end{eqnarray*}

""" sinhcosh!


doc"""
    tanh!(c, a, p, k)

Update the `k-th` expansion coefficients `c[k+1]` of `c = tanh(a)`,
for `c` and `a` either `Taylor1` or `TaylorN`; `p = a^2` and
is passed as an argument for efficiency.

\begin{equation*}
c_k = a_k - \frac{1}{k} \sum_{j=0}^{k-1} (k-j) a_{k-j} p_j.
\end{equation*}

""" tanh!
