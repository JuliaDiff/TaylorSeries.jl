# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

## Broadcast for Taylor1 and TaylorN

import .Broadcast: BroadcastStyle, Broadcasted, broadcasted

# BroadcastStyle definitions and basic precedence rules
struct Taylor1Style{T} <: Base.Broadcast.AbstractArrayStyle{0} end
Taylor1Style{T}(::Val{N}) where {T, N}= Base.Broadcast.DefaultArrayStyle{N}()
BroadcastStyle(::Type{<:Taylor1{T}}) where {T} = Taylor1Style{T}()
BroadcastStyle(::Taylor1Style{T}, ::Base.Broadcast.DefaultArrayStyle{0}) where {T} = Taylor1Style{T}()
BroadcastStyle(::Taylor1Style{T}, ::Base.Broadcast.DefaultArrayStyle{1}) where {T} = Base.Broadcast.DefaultArrayStyle{1}()
#
struct HomogeneousPolynomialStyle{T} <: Base.Broadcast.AbstractArrayStyle{0} end
HomogeneousPolynomialStyle{T}(::Val{N}) where {T, N}= Base.Broadcast.DefaultArrayStyle{N}()
BroadcastStyle(::Type{<:HomogeneousPolynomial{T}}) where {T} = HomogeneousPolynomialStyle{T}()
BroadcastStyle(::HomogeneousPolynomialStyle{T}, ::Base.Broadcast.DefaultArrayStyle{0}) where {T} = HomogeneousPolynomialStyle{T}()
BroadcastStyle(::HomogeneousPolynomialStyle{T}, ::Base.Broadcast.DefaultArrayStyle{1}) where {T} = Base.Broadcast.DefaultArrayStyle{1}()
#
struct TaylorNStyle{T} <: Base.Broadcast.AbstractArrayStyle{0} end
TaylorNStyle{T}(::Val{N}) where {T, N}= Base.Broadcast.DefaultArrayStyle{N}()
BroadcastStyle(::Type{<:TaylorN{T}}) where {T} = TaylorNStyle{T}()
BroadcastStyle(::TaylorNStyle{T}, ::Base.Broadcast.DefaultArrayStyle{0}) where {T} = TaylorNStyle{T}()
BroadcastStyle(::TaylorNStyle{T}, ::Base.Broadcast.DefaultArrayStyle{1}) where {T} = Base.Broadcast.DefaultArrayStyle{1}()


# Extend eltypes so things like [1.0] .+ t work
Base.Broadcast.eltypes(t::Tuple{Taylor1,AbstractArray}) =
    Tuple{Base.Broadcast._broadcast_getindex_eltype([t[1]]), Base.Broadcast._broadcast_getindex_eltype(t[2])}
Base.Broadcast.eltypes(t::Tuple{AbstractArray,Taylor1}) =
    Tuple{Base.Broadcast._broadcast_getindex_eltype(t[1]), Base.Broadcast._broadcast_getindex_eltype([t[2]])}
Base.Broadcast.eltypes(t::Tuple{HomogeneousPolynomial,AbstractArray}) =
    Tuple{Base.Broadcast._broadcast_getindex_eltype([t[1]]), Base.Broadcast._broadcast_getindex_eltype(t[2])}
Base.Broadcast.eltypes(t::Tuple{AbstractArray,HomogeneousPolynomial}) =
    Tuple{Base.Broadcast._broadcast_getindex_eltype(t[1]), Base.Broadcast._broadcast_getindex_eltype([t[2]])}
Base.Broadcast.eltypes(t::Tuple{TaylorN,AbstractArray}) =
    Tuple{Base.Broadcast._broadcast_getindex_eltype([t[1]]), Base.Broadcast._broadcast_getindex_eltype(t[2])}
Base.Broadcast.eltypes(t::Tuple{AbstractArray,TaylorN}) =
    Tuple{Base.Broadcast._broadcast_getindex_eltype(t[1]), Base.Broadcast._broadcast_getindex_eltype([t[2]])}

# We follow https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-iteration-1
"`A = find_taylor(As)` returns the first Taylor1 among the arguments."
find_taylor(bc::Broadcasted) = find_taylor(bc.args)
find_taylor(args::Tuple) = find_taylor(find_taylor(args[1]), Base.tail(args))
find_taylor(x) = x
find_taylor(a::Taylor1, rest) = a
find_taylor(a::HomogeneousPolynomial, rest) = a
find_taylor(a::TaylorN, rest) = a
find_taylor(::AbstractArray, rest) = find_taylor(rest)

# Extend Base.similar
function Base.similar(bc::Broadcasted{Taylor1Style{S}}, ::Type{T}) where {S, T}
    # Proper promotion
    R = Base.Broadcast.combine_eltypes(bc.f, bc.args)
    # Scan the inputs for the Taylor1:
    A = find_taylor(bc)
    # Create the output
    return Taylor1(similar(A.coeffs, R), A.order)
end

function Base.similar(bc::Broadcasted{HomogeneousPolynomialStyle{S}}, ::Type{T}) where {S, T}
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

function Base.similar(bc::Broadcasted{TaylorNStyle{S}}, ::Type{T}) where {S, T}
    # Proper promotion
    R = Base.Broadcast.combine_eltypes(bc.f, bc.args)
    # Scan the inputs for the TaylorN:
    A = find_taylor(bc)
    # Create the output
    return TaylorN(similar(A.coeffs, R), A.order)
end


# Adapted from Base.Broadcast.copyto!, base/broadcasting.jl, line 832
@inline function Base.copyto!(dest::Taylor1{T}, bc::Broadcasted) where {T}
    axes(dest) == axes(bc) || Base.Broadcast.throwdm(axes(dest), axes(bc))
    # Performance optimization: broadcast!(identity, dest, A) is equivalent to copyto!(dest, A) if indices match
    if bc.f === identity && bc.args isa Tuple{Taylor1{T}} # only a single input argument to broadcast!
        A = bc.args[1]
        if axes(dest) == axes(A)
            return copyto!(dest, A)
        end
    end
    bc′ = Base.Broadcast.preprocess(dest, bc)
    copyto!(dest, bc′[1])
    return dest
end

# @inline function Base.copyto!(dest::HomogeneousPolynomial{T}, bc::Broadcasted) where {T<:NumberNotSeries}
#     axes(dest) == axes(bc) || Base.Broadcast.throwdm(axes(dest), axes(bc))
#     # Performance optimization: broadcast!(identity, dest, A) is equivalent to copyto!(dest, A) if indices match
#     if bc.f === identity && bc.args isa Tuple{AbstractArray} # only a single input argument to broadcast!
#         A = bc.args[1]
#         if axes(dest) == axes(A)
#             return copyto!(dest, A)
#         end
#     end
#     bc′= Base.Broadcast.preprocess(dest.coeffs, bc)
#     # I is the coefficients index
#     @simd for I in eachindex(bc′)
#         # @inbounds dest[I] = getcoeff(bc′[I], I)
#         @inbounds dest[I] = bc′[I].coeffs[I]
#     end
#     return dest
# end

@inline function Base.copyto!(dest::TaylorN{T}, bc::Broadcasted) where {T<:NumberNotSeries}
    axes(dest) == axes(bc) || Base.Broadcast.throwdm(axes(dest), axes(bc))
    # Performance optimization: broadcast!(identity, dest, A) is equivalent to copyto!(dest, A) if indices match
    if bc.f === identity && bc.args isa Tuple{TaylorN{T}} # only a single input argument to broadcast!
        A = bc.args[1]
        if axes(dest) == axes(A)
            return copyto!(dest, A)
        end
    end
    bc′ = Base.Broadcast.preprocess(dest.coeffs, bc)
    # I is the coefficients index
    @simd for I in eachindex(bc′.args)
        @show(I, bc′.args[I], dest[I])
        @inbounds dest[I] = getcoeff(bc′[I], I)
    end
    return dest
end


# Broadcasted extensions
# This prevents broadcasting being applied to `a` and `b`
# for the mutating functions, and to act only in `k`
function broadcasted(::Taylor1Style{T}, f!, a::Taylor1{T}, b::Taylor1{T}, k) where {T}
    @inbounds for i in eachindex(k)
        f!(a, b, k[i])
    end
    nothing
end
