using .IntervalArithmetic

"""
    evaluate(a, dx; [horner], [normalization], [intersection])

Evaluate a `Taylor1` polynomial with interval coefficients.

Use the `horner` argument (`true` by default), to apply Horner's rule.

If `normalization` is different from ("none" the default), evaluate after having
re-scaled `a` and `dx` to: (i) the interval `[-1, 1]` if `normalization="symmetric"`,
(ii) the unit interval `[0, 1]` if `normalization="asymmetric"`.

To combine the previous arguments, use the `intersection` option (`false` by default).
This method applies the four possible combinations of `horner` and `normalize`, and
returns their intersection computed with interval arithmetic.
In this way, the tightest possible bound with any of the available methods is obtained.
"""
function evaluate(a::Taylor1{<:Interval{T}}, dx::Interval{T};
                  horner=true, normalization="none", intersection=false) where {T<:Number}
    # if intersection is true, *ignore* other options
    if intersection
        I1 = evaluate(a, dx, horner=true, normalization="none")
        I2 = evaluate(a, dx, horner=false, normalization="none")
        I3 = evaluate(a, dx, horner=true, normalization="symmetric")
        I4 = evaluate(a, dx, horner=false, normalization="symmetric")
        I5 = evaluate(a, dx, horner=true, normalization="asymmetric")
        I6 = evaluate(a, dx, horner=false, normalization="asymmetric")
        return intersect(I1, I2, I3, I4, I5, I6)
    end

    if normalization != "none"
        symmetric = normalization == "symmetric" ? true : false
        (a, dx) = _normalize(a, dx; symmetric=symmetric)
    end

    if horner
        suma = _evaluate_horner(a, dx)
    else
        suma = _evaluate_naive(a, dx)
    end
    return suma
end

function _evaluate_horner(a, dx)
    @inbounds begin
        suma = a[end]
        for k in a.order-1:-1:0
            suma = suma*dx + a[k]
        end
    end
    return suma
end

function _evaluate_naive(a, dx)
    @inbounds begin
        suma = a[0]
        for k in 1:a.order
            suma = suma + a[k] * dx^k
        end
    end
    return suma
end


"""
    normalize_interval(a::Taylor1, I::Interval, symI::Bool=true)

Normalizes `a::Taylor1` such that the interval `I` is mapped
by an affine transformation to the interval `-1..1` (`symI=true`)
or to `0..1` (`symI=false`).
"""
normalize_taylor(a::Taylor1, I::Interval{T}, symI::Bool=true) where {T} =
    _normalize(a, I, Val(symI))

"""
    normalize_interval(a::TaylorN, I::IntervalBox, symI::Bool=true)

Normalize `a::TaylorN` such that the intervals in `I::IntervalBox`
are mapped by an affine transformation to the intervals `-1..1`
(`symI=true`) or to `0..1` (`symI=false`).
"""
normalize_taylor(a::TaylorN, I::IntervalBox{N,T}, symI::Bool=true) where {N,T} =
    _normalize(a, I, Val(symI))


#  I -> -1..1)
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
    x = Array{typeof(a)}(undef, N)
    for ind in eachindex(x)
        x[ind] = TaylorN(ind, order=order)
    end
    x = mid.(I) .+ x .* radius.(I)
    return a(x)
end
#  I -> IntervalBox(0..1, Val(N))
function _normalize(a::TaylorN, I::IntervalBox{N,T}, ::Val{false}) where {N,T}
    order = get_order(a)
    x = Array{typeof(a),1}(undef, N)
    for ind in eachindex(x)
        x[ind] = TaylorN(ind, order=order)
    end
    x .= inf.(I) .+ x .* diam.(I)
    return a(x)
end
