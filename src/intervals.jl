using .IntervalArithmetic

# """
#     evaluate(a, dx; [horner], [normalization], [intersection])
#
# Evaluate a `Taylor1` polynomial with interval coefficients.
#
# Use the `horner` argument (`true` by default), to apply Horner's rule.
#
# If `normalization` is different from ("none" the default), evaluate after having
# re-scaled `a` and `dx` to: (i) the interval `[-1, 1]` if `normalization="symmetric"`,
# (ii) the unit interval `[0, 1]` if `normalization="asymmetric"`.
#
# To combine the previous arguments, use the `intersection` option (`false` by default).
# This method applies the four possible combinations of `horner` and `normalize`, and
# returns their intersection computed with interval arithmetic.
# In this way, the tightest possible bound with any of the available methods is obtained.
# """
# function evaluate(a::Taylor1{<:Interval{T}}, dx::Interval{T};
#                   horner=true, normalization="none", intersection=false) where {T<:Number}
#     # if intersection is true, *ignore* other options
#     if intersection
#         I1 = evaluate(a, dx, horner=true, normalization="none")
#         I2 = evaluate(a, dx, horner=false, normalization="none")
#         I3 = evaluate(a, dx, horner=true, normalization="symmetric")
#         I4 = evaluate(a, dx, horner=false, normalization="symmetric")
#         I5 = evaluate(a, dx, horner=true, normalization="asymmetric")
#         I6 = evaluate(a, dx, horner=false, normalization="asymmetric")
#         return intersect(I1, I2, I3, I4, I5, I6)
#     end
#
#     if normalization != "none"
#         symmetric = normalization == "symmetric" ? true : false
#         (a, dx) = _normalize(a, dx; symmetric=symmetric)
#     end
#
#     if horner
#         suma = _evaluate_horner(a, dx)
#     else
#         suma = _evaluate_naive(a, dx)
#     end
#     return suma
# end
#
# function _evaluate_horner(a, dx)
#     @inbounds begin
#         suma = a[end]
#         for k in a.order-1:-1:0
#             suma = suma*dx + a[k]
#         end
#     end
#     return suma
# end
#
# function _evaluate_naive(a, dx)
#     @inbounds begin
#         suma = a[0]
#         for k in 1:a.order
#             suma = suma + a[k] * dx^k
#         end
#     end
#     return suma
# end


function evaluate(a::Taylor1, dx::Interval)
    dx == (-1..1) && return _evaluate(a, dx, Val(true))
    dx == (0..1) && return _evaluate(a, dx, Val(false))
    # Usual Horner rule
    @inbounds begin
        suma = a[end]*one(dx)
        for k in a.order-1:-1:0
            suma = suma*dx + a[k]
        end
    end
    return suma
end

function _evaluate(a::Taylor1, dx::Interval, ::Val{true})
    if iseven(a.order)
        @inbounds begin
            suma1 = a[end-1]
            for k in a.order-1:-2:1
                suma1 += a[k]
            end
            suma2 = a[end]
            for k in a.order-2:-2:1
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
            for k in a.order-1:-2:1
                suma2 += a[k]
            end
        end
    end
    return a[0] + suma1*dx + suma2 * (one(dx) * (0..1))
end

function _evaluate(a::Taylor1, dx::Interval, ::Val{false})
    @inbounds begin
        suma = a[end]
        for k in a.order-1:-1:1
            suma = suma + a[k]
        end
    end
    return a[0] + suma*dx
end


# function evaluate(a::TaylorN, dx::IntervalBox{N,T}) where {N,T}
#
#     @assert N == get_numvars()
#
#     @show(dx == IntervalBox(-1..1, Val(N))) && return _evaluate(a, dx, Val(true))
#     @show(dx == IntervalBox(0..1, Val(N))) && return _evaluate(a, dx, Val(false))
#
#     # Otherwise, usual method
#     R = promote_type(T,S)
#     a_length = length(a)
#     suma = zeros(R, a_length)
#     for homPol in 1:length(a)
#             # sun = zero(R)
#             # for (i, a_coeff) in enumerate(a.coeffs[homPol].coeffs)
#             #     iszero(a_coeff) && continue
#             #     tmp = vals[1]^(coeff_table[homPol][i][1])
#             #     for n in 2:N
#             #         tmp *= vals[n]^(coeff_table[homPol][i][n])
#             #     end
#             #     sun += a_coeff * tmp
#             # end
#             suma[homPol] = a.coeffs[homPol](vals)
#         end
#
#         return sum( sort!(suma, by=abs2) )
#     end
# end
#
# function _evaluate(a::TaylorN, dx::IntervalBox{N,T}, ::Val{true}) where {N,T}
#     dx_even = one(dx) * IntervalBox(0..1, Val(N))
#     if iseven(a.order)
#         @inbounds begin
#             suma = a[end]*dx_even
#             for k in a.order-1:-2:1
#                 suma = suma + a[k] * dx
#             end
#             for k in a.order-2:-2:0
#                 suma = suma + a[k] * dx_even
#             end
#         end
#     else
#         @inbounds begin
#             suma = a[end]*dx
#             for k in a.order-2:-2:1
#                 suma = suma + a[k] * dx
#             end
#             for k in a.order-1:-2:0
#                 suma = suma + a[k] * dx_even
#             end
#         end
#     end
#     return suma
# end
#
# function _evaluate(a::Taylor1, dx::Interval, ::Val{false})
#     @inbounds begin
#         suma = a[end]*dx
#         for k in a.order-1:-1:0
#             suma = suma + a[k]*dx
#         end
#     end
#     return suma
# end

function evaluate(a::TaylorN, dx::IntervalBox{N,T}) where {N,T}

    @assert N == get_numvars()
    a_length = length(a)
    suma = constant_term(a) + Interval{T}(0, 0)
    for homPol in 1:length(a)
        suma += evaluate(a.coeffs[homPol], dx)
    end

    return suma
end

function evaluate(a::HomogeneousPolynomial, dx::IntervalBox{N,T}) where {T, N}

    @assert N == get_numvars()
    dx == IntervalBox(-1..1, Val(N)) && return _evaluate(a, dx, Val(true))
    dx == IntervalBox(0..1, Val(N)) && return _evaluate(a, dx, Val(false))

    return evaluate(a, dx...)
end

function _evaluate(a::HomogeneousPolynomial, dx::IntervalBox{N,T}, ::Val{true} ) where {T, N}
    a.order == 0 && return a[1] + Interval{T}(0, 0)

    suma = zero(a[1])
    for homPol in a.coeffs
        suma += homPol
    end
    iseven(a.order) && return suma * Interval{T}(0, 1)
    return suma * dx[1]
end

function _evaluate(a::HomogeneousPolynomial, dx::IntervalBox{N,T}, ::Val{false} ) where {T, N}
    a.order == 0 && return a[1] + Interval{T}(0, 0)

    suma = zero(a[1])
    for homPol in a.coeffs
        suma += homPol
    end
    return suma * dx[1]
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
