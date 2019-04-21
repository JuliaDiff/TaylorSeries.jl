# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

# Broadcasting: Implemented for Taylor1 and TaylorN
import .Broadcast: BroadcastStyle, Style, AbstractArrayStyle, Broadcasted, broadcasted

# BroadcastStyle definitions
BroadcastStyle(::Type{<:Taylor1{T}}) where {T}= Style{Taylor1{T}}()
BroadcastStyle(::Style{Taylor1{T}}, ::AbstractArrayStyle{0}) where {T}= Style{Taylor1{T}}()
BroadcastStyle(::Style{Taylor1{T}}, ::AbstractArrayStyle{1}) where {T}= Style{Taylor1{T}}()
#
BroadcastStyle(::Type{<:HomogeneousPolynomial{T}}) where {T} = Style{HomogeneousPolynomial{T}}()
BroadcastStyle(::Style{HomogeneousPolynomial{T}}, ::AbstractArrayStyle{0}) where {T} = Style{HomogeneousPolynomial{T}}()
BroadcastStyle(::Style{HomogeneousPolynomial{T}}, ::AbstractArrayStyle{1}) where {T} = Style{HomogeneousPolynomial{T}}()
#
BroadcastStyle(::Type{<:TaylorN{T}}) where {T} = Style{TaylorN{T}}()
BroadcastStyle(::Style{TaylorN{T}}, ::AbstractArrayStyle{0}) where {T} = Style{TaylorN{T}}()
BroadcastStyle(::Style{TaylorN{T}}, ::AbstractArrayStyle{1}) where {T} = Style{TaylorN{T}}()

# Precedence rules (for mixtures)
BroadcastStyle(::Style{Taylor1{Taylor1{T}}}, ::Style{Taylor1{S}}) where
    {T,S} = Style{Taylor1{Taylor1{T}}}()
#
BroadcastStyle(::Style{Taylor1{TaylorN{T}}}, ::Style{HomogeneousPolynomial{S}}) where
    {T<:NumberNotSeries, S<:NumberNotSeries} = Style{Taylor1{TaylorN{promote_type(T,S)}}}()
BroadcastStyle(::Style{Taylor1{TaylorN{T}}}, ::Style{TaylorN{S}}) where
    {T<:NumberNotSeries, S<:NumberNotSeries} = Style{Taylor1{TaylorN{promote_type(T,S)}}}()
BroadcastStyle(::Style{Taylor1{T}}, ::Style{HomogeneousPolynomial{Taylor1{S}}}) where
    {T<:NumberNotSeries, S<:NumberNotSeries} = Style{HomogeneousPolynomial{Taylor1{promote_type(T,S)}}}()
BroadcastStyle(::Style{Taylor1{T}}, ::Style{TaylorN{Taylor1{S}}}) where
    {T<:NumberNotSeries, S<:NumberNotSeries} = Style{TaylorN{Taylor1{promote_type(T,S)}}}()


# We follow https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-iteration-1
"`A = find_taylor(As)` returns the first Taylor1 among the arguments."
find_taylor(bc::Broadcasted) = find_taylor(bc.args)
find_taylor(args::Tuple) = find_taylor(find_taylor(args[1]), Base.tail(args))
find_taylor(x) = x
find_taylor(a::Taylor1, rest) = a
find_taylor(a::HomogeneousPolynomial, rest) = a
find_taylor(a::TaylorN, rest) = a
find_taylor(::Any, rest) = find_taylor(rest)

# Extend Base.similar
function Base.similar(bc::Broadcasted{Style{Taylor1{S}}}, ::Type{T}) where {S, T}
    # Proper promotion
    R = Base.Broadcast.combine_eltypes(bc.f, bc.args)
    # Scan the inputs for the Taylor1:
    A = find_taylor(bc)
    # Create the output
    return Taylor1(similar(A.coeffs, R), A.order)
end

function Base.similar(bc::Broadcasted{Style{HomogeneousPolynomial{S}}}, ::Type{T}) where {S, T}
    # Proper promotion
    @show(Base.Broadcast.eltypes(bc.args))
    # combine_eltypes(f, args::Tuple) = Base._return_type(f, eltypes(args))
    R = Base.Broadcast.combine_eltypes(bc.f, bc.args)
    @show(bc.f, bc.args)
    # Scan the inputs for the HomogeneousPolynomial:
    A = find_taylor(bc)
    @show(A, R)
    # Create the output
    return HomogeneousPolynomial(similar(A.coeffs, R), A.order)
end

function Base.similar(bc::Broadcasted{Style{TaylorN{S}}}, ::Type{T}) where {S, T}
    # Proper promotion
    R = Base.Broadcast.combine_eltypes(bc.f, bc.args)
    # Scan the inputs for the TaylorN:
    A = find_taylor(bc)
    # Create the output
    return TaylorN(similar(A.coeffs, R), A.order)
end


# Adapted from Base.Broadcast.copyto!, base/broadcasting.jl, line 832
@inline function Base.copyto!(dest::Taylor1{T}, bc::Broadcasted) where {T<:NumberNotSeries}
    axes(dest) == axes(bc) || Base.Broadcast.throwdm(axes(dest), axes(bc))
    # Performance optimization: broadcast!(identity, dest, A) is equivalent to copyto!(dest, A) if indices match
    if bc.f === identity && bc.args isa Tuple{AbstractArray} # only a single input argument to broadcast!
        A = bc.args[1]
        if axes(dest) == axes(A)
            return copyto!(dest, A)
        end
    end
    bc′ = Base.Broadcast.preprocess(dest.coeffs, bc)
    # I is the coefficients index
    @simd for I in eachindex(bc′)
        @inbounds dest[I-1] = getcoeff(bc′[I], I-1)
    end
    return dest
end

@inline function Base.copyto!(dest::HomogeneousPolynomial{T}, bc::Broadcasted) where {T<:NumberNotSeries}
    @show(dest, get_order(dest))
    axes(dest) == axes(bc) || Base.Broadcast.throwdm(axes(dest), axes(bc))
    # Performance optimization: broadcast!(identity, dest, A) is equivalent to copyto!(dest, A) if indices match
    if bc.f === identity && bc.args isa Tuple{AbstractArray} # only a single input argument to broadcast!
        A = bc.args[1]
        if axes(dest) == axes(A)
            return copyto!(dest, A)
        end
    end
    bc′= Base.Broadcast.preprocess(dest.coeffs, bc)
    @show(bc′)
    # I is the coefficients index
    @simd for I in eachindex(bc′)
        # @show(I, bc′[I].coeffs, bc′[I].order)
        # @inbounds dest[I] = getcoeff(bc′[I], I)
        @inbounds dest[I] = bc′[I].coeffs[I]
    end
    return dest
end

@inline function Base.copyto!(dest::TaylorN{T}, bc::Broadcasted) where {T<:NumberNotSeries}
    @show(dest, bc)
    axes(dest) == axes(bc) || Base.Broadcast.throwdm(axes(dest), axes(bc))
    # Performance optimization: broadcast!(identity, dest, A) is equivalent to copyto!(dest, A) if indices match
    if bc.f === identity && bc.args isa Tuple{AbstractArray} # only a single input argument to broadcast!
        A = bc.args[1]
        if axes(dest) == axes(A)
            return copyto!(dest, A)
        end
    end
    bc′ = Base.Broadcast.preprocess(dest.coeffs, bc)
    # I is the coefficients index
    @simd for I in eachindex(bc′)
        @show(I, bc´[I])
        @inbounds dest[I] = getcoeff(bc′[I], I)
    end
    return dest
end


# Broadcasted extensions
# This should prevent broadcasting being applied in `a` and `b`
# for the mutating functions, and work only in `k`
function broadcasted(::Style{Taylor1{T}}, f!, a::Taylor1{T}, b::Taylor1{T}, k) where {T}
    @inbounds for i in eachindex(k)
        f!(a, b, k[i])
    end
    nothing
end
