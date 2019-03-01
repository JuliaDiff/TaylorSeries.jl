using .IntervalArithmetic

function evaluate(a::Taylor1, dx::Interval)
    dx == (-1..1) && return _evaluate(a, dx, Val(true))
    dx == (0..1) && return _evaluate(a, dx, Val(false))
    # Usual Horner rule
    uno = one(dx)
    @inbounds begin
        suma = a[end]*uno
        for k in a.order-1:-1:0
            suma = suma*dx + a[k]*uno
        end
    end
    return suma
end

function _evaluate(a::Taylor1, dx::Interval{T}, ::Val{true}) where {T}
    if iseven(a.order)
        @inbounds begin
            suma1 = a[end-1]
            for k in a.order-3:-2:1
                suma1 += a[k]
            end
            suma2 = a[end]
            for k in a.order-2:-2:2
                suma2 += a[k]
            end
        end
    else
        @inbounds begin
            suma1 = a[end]
            for k in a.order-2:-2:1
                suma1 += a[k]
            end
            suma2 = a[end-1]
            for k in a.order-3:-2:2
                suma2 += a[k]
            end
        end
    end
    return a[0] + suma1*dx + suma2 * Interval{T}(0, 1)
end

function _evaluate(a::Taylor1, dx::Interval, ::Val{false})
    uno = one(dx)
    @inbounds begin
        suma = a[end]*uno
        for k in a.order-1:-1:1
            suma = suma*dx + a[k]*uno
        end
    end
    return a[0] + suma
end


function evaluate(a::TaylorN, dx::IntervalBox{N,T}) where {N,T}

    @assert N == get_numvars()
    a_length = length(a)
    suma = zero(constant_term(a)) + Interval{T}(0, 0)
    @inbounds for homPol in length(a):-1:1
        suma += evaluate(a.coeffs[homPol], dx)
    end

    return suma
end

function evaluate(a::HomogeneousPolynomial, dx::IntervalBox{N,T}) where {T, N}

    @assert N == get_numvars()
    dx == IntervalBox(-1..1, Val(N)) && return _evaluate(a, dx, Val(true))
    dx == IntervalBox( 0..1, Val(N)) && return _evaluate(a, dx, Val(false))

    return evaluate(a, dx...)
end

function _evaluate(a::HomogeneousPolynomial, dx::IntervalBox{N,T}, ::Val{true} ) where {T, N}
    a.order == 0 && return a[1] + Interval{T}(0, 0)

    ct = coeff_table[a.order+1]
    @inbounds suma = a[1]*Interval{T}(0,0)

    Ieven = Interval{T}(0,1)
    for (i,a_coeff) in enumerate(a.coeffs)
        iszero(a_coeff) && continue
        if isodd(sum(ct[i]))
            suma += sum(a_coeff) * dx[1]
            continue
        end
        @inbounds tmp = iseven(ct[i][1]) ? Ieven : dx[1]
        for n in 2:N
            @inbounds vv = iseven(ct[i][n]) ? Ieven : dx[1]
            tmp *= vv
        end
        suma += a_coeff * tmp
    end
    return suma
end

function _evaluate(a::HomogeneousPolynomial, dx::IntervalBox{N,T}, ::Val{false} ) where {T, N}
    a.order == 0 && return a[1] + Interval{T}(0, 0)

    @inbounds suma = zero(a[1])*dx[1]
    @inbounds for homPol in a.coeffs
        suma += homPol*dx[1]
    end
    return suma
end


"""
    normalize_taylor(a::Taylor1, I::Interval, symI::Bool=true)

Normalizes `a::Taylor1` such that the interval `I` is mapped
by an affine transformation to the interval `-1..1` (`symI=true`)
or to `0..1` (`symI=false`).
"""
normalize_taylor(a::Taylor1, I::Interval{T}, symI::Bool=true) where {T} =
    _normalize(a, I, Val(symI))

"""
    normalize_taylor(a::TaylorN, I::IntervalBox, symI::Bool=true)

Normalize `a::TaylorN` such that the intervals in `I::IntervalBox`
are mapped by an affine transformation to the intervals `-1..1`
(`symI=true`) or to `0..1` (`symI=false`).
"""
normalize_taylor(a::TaylorN, I::IntervalBox{N,T}, symI::Bool=true) where {N,T} =
    _normalize(a, I, Val(symI))


#  I -> -1..1
function _normalize(a::Taylor1, I::Interval{T}, ::Val{true}) where {T}
    order = get_order(a)
    t = Taylor1(T, order)
    tnew = mid(I) + t*radius(I)
    return a(tnew)
end

#  I -> 0..1
function _normalize(a::Taylor1, I::Interval{T}, ::Val{false}) where {T}
    order = get_order(a)
    t = Taylor1(T, order)
    tnew = inf(I) + t*diam(I)
    return a(tnew)
end

#  I -> IntervalBox(-1..1, Val(N))
function _normalize(a::TaylorN, I::IntervalBox{N,T}, ::Val{true}) where {N,T}
    order = get_order(a)
    x = Vector{typeof(a)}(undef, N)
    for ind in eachindex(x)
        x[ind] = mid(I[ind]) + TaylorN(ind, order=order)*radius(I[ind])
    end
    return a(x)
end

#  I -> IntervalBox(0..1, Val(N))
function _normalize(a::TaylorN, I::IntervalBox{N,T}, ::Val{false}) where {N,T}
    order = get_order(a)
    x = Vector{typeof(a)}(undef, N)
    for ind in eachindex(x)
        x[ind] = inf(I[ind]) + TaylorN(ind, order=order)*diam(I[ind])
    end
    return a(x)
end
