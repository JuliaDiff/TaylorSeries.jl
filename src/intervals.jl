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

function _normalize(a::Taylor1{<:Interval{T}}, dx::Interval{T}; symmetric=false) where {T<:Number}
    t = Taylor1(T, a.order)
    if symmetric
        tnew = t*diam(dx)/2 + mid(dx)
        Inew = Interval{T}(-1, 1)
    else
        tnew = t*diam(dx) + dx.lo
        Inew = Interval{T}(0, 1)
    end
    return (a(tnew), Inew)
end
