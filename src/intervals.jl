using .IntervalArithmetic

"""
    evaluate(a, dx; [horner], [normalize], [intersect])

Evaluate a `Taylor1` polynomial with interval coefficients.

Use the `horner` argument (`true` by default), to apply Horner's rule.

Use the `normalize` argument (`false` by default) to evaluate after having
re-scaled the inputs to the unit symmetric interval `[-1, 1]`.

To combine the previous arguments, use the `intersect` option (`false` by default).
This method applies the four possible combinations of `horner` and `normalize`, and
returns their intersection computed with interval arithmetic.
In this way, the tightest possible bound with any of the available methods is obtained.
"""
function evaluate(a::Taylor1{<:Interval{T}}, dx::Interval{T};
                  horner=true, normalize=false, intersect=false) where {T<:Number}
    # if intersect is true, *ignore* other options
    if intersect
        return evaluate(a, dx, horner=false, normalize=false) ∩
               evaluate(a, dx, horner=false, normalize=true) ∩
               evaluate(a, dx, horner=true, normalize=false) ∩
               evaluate(a, dx, horner=true, normalize=true)
    end

    if normalize
        (a, dx) = _normalize(a, dx)
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

function _normalize(a::Taylor1{<:Interval{T}}, dx::Interval{T}) where {T<:Number}
    t = (0..0) + Taylor1(T, a.order)
    return (a(radius(dx) + mid(dx)*t), Interval{T}(-1, 1))
end
