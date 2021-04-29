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
Taylor1Style{T}(::Val{N}) where {T,N}= Base.Broadcast.DefaultArrayStyle{N}()
BroadcastStyle(::Type{<:Taylor1{T}}) where {T} = Taylor1Style{T}()
BroadcastStyle(::Taylor1Style{T}, ::Base.Broadcast.DefaultArrayStyle{0}) where {T} = Taylor1Style{T}()
BroadcastStyle(::Taylor1Style{T}, ::Base.Broadcast.DefaultArrayStyle{1}) where {T} = Base.Broadcast.DefaultArrayStyle{1}()
#
struct HomogeneousPolynomialStyle{T} <: Base.Broadcast.AbstractArrayStyle{0} end
HomogeneousPolynomialStyle{T}(::Val{N}) where {T,N}= Base.Broadcast.DefaultArrayStyle{N}()
BroadcastStyle(::Type{<:HomogeneousPolynomial{T}}) where {T} = HomogeneousPolynomialStyle{T}()
BroadcastStyle(::HomogeneousPolynomialStyle{T}, ::Base.Broadcast.DefaultArrayStyle{0}) where {T} = HomogeneousPolynomialStyle{T}()
BroadcastStyle(::HomogeneousPolynomialStyle{T}, ::Base.Broadcast.DefaultArrayStyle{1}) where {T} = Base.Broadcast.DefaultArrayStyle{1}()
#
struct TaylorNStyle{T} <: Base.Broadcast.AbstractArrayStyle{0} end
TaylorNStyle{T}(::Val{N}) where {T, N}= Base.Broadcast.DefaultArrayStyle{N}()
BroadcastStyle(::Type{<:TaylorN{T}}) where {T} = TaylorNStyle{T}()
BroadcastStyle(::TaylorNStyle{T}, ::Base.Broadcast.DefaultArrayStyle{0}) where {T} = TaylorNStyle{T}()
BroadcastStyle(::TaylorNStyle{T}, ::Base.Broadcast.DefaultArrayStyle{1}) where {T} = Base.Broadcast.DefaultArrayStyle{1}()

# Precedence rules for mixtures
BroadcastStyle(::TaylorNStyle{Taylor1{T}}, ::Taylor1Style{T}) where {T} = TaylorNStyle{Taylor1{T}}()
BroadcastStyle(::Taylor1Style{TaylorN{T}}, ::TaylorNStyle{T}) where {T} = Taylor1Style{TaylorN{T}}()

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


# # We follow https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-iteration-1
# "`A = find_taylor(As)` returns the first Taylor1 among the arguments."
# find_taylor(bc::Broadcasted) = find_taylor(bc.args)
# find_taylor(args::Tuple) = find_taylor(find_taylor(args[1]), Base.tail(args))
# find_taylor(x) = x
# find_taylor(a::Taylor1, rest) = a
# find_taylor(a::HomogeneousPolynomial, rest) = a
# find_taylor(a::TaylorN, rest) = a
# find_taylor(::AbstractArray, rest) = find_taylor(rest)
#
# # Extend similar
# function similar(bc::Broadcasted{Taylor1Style{S}}, ::Type{T}) where {S, T}
#     # Proper promotion
#     R = Base.Broadcast.combine_eltypes(bc.f, bc.args)
#     # Scan the inputs for the Taylor1:
#     A = find_taylor(bc)
#     # Create the output
#     return Taylor1(similar(A.coeffs, R), A.order)
# end
#
# function similar(bc::Broadcasted{HomogeneousPolynomialStyle{S}}, ::Type{T}) where {S, T}
#     # Proper promotion
#     # combine_eltypes(f, args::Tuple) = Base._return_type(f, eltypes(args))
#     R = Base.Broadcast.combine_eltypes(bc.f, bc.args)
#     # Scan the inputs for the HomogeneousPolynomial:
#     A = find_taylor(bc)
#     # Create the output
#     return HomogeneousPolynomial(similar(A.coeffs, R), A.order)
# end
#
# function similar(bc::Broadcasted{TaylorNStyle{S}}, ::Type{T}) where {S, T}
#     # Proper promotion
#     R = Base.Broadcast.combine_eltypes(bc.f, bc.args)
#     # Scan the inputs for the TaylorN:
#     A = find_taylor(bc)
#     # Create the output
#     return TaylorN(similar(A.coeffs, R), A.order)
# end


# Adapted from Base.Broadcast.copyto!, base/broadcasting.jl, line 832
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval begin
        @inline function copyto!(dest::$T{T}, bc::Broadcasted) where {T<:Number}
            axes(dest) == axes(bc) || Base.Broadcast.throwdm(axes(dest), axes(bc))
            # Performance optimization: broadcast!(identity, dest, A) is equivalent to copyto!(dest, A) if indices match
            if bc.f === identity && bc.args isa Tuple{$T{T}} # only a single input argument to broadcast!
                A = bc.args[1]
                if axes(dest) == axes(A)
                    return copyto!(dest, A)
                end
            end
            bc′ = Base.Broadcast.preprocess(dest, bc)
            copyto!(dest, bc′[1])
            return dest
        end
    end
end


# Broadcasted extensions
@inline broadcasted(::Taylor1Style{T}, ::Type{Float32}, a::Taylor1{T}) where {T<:Number} =
    Taylor1(Float32.(a.coeffs), a.order)
@inline broadcasted(::TaylorNStyle{T}, ::Type{Float32}, a::TaylorN{T}) where {T<:Number} =
    convert(TaylorN{Float32}, a)

# # This prevents broadcasting being applied to the Taylor1/TaylorN params
# # for the mutating functions, and to act only in `k`
# for (T, TS) in ((:Taylor1, :Taylor1Style), (:TaylorN, :TaylorNStyle))
#     for f in (add!, subst!, sqr!, sqrt!, exp!, log!, identity!, zero!,
#         one!, abs!, abs2!, deg2rad!, rad2deg!)
#         @eval begin
#             @inline function broadcasted(::$TS{T}, fn::typeof($f), r::$T{T}, a::$T{T}, k) where {T}
#                 @inbounds for i in eachindex(k)
#                     fn(r, a, k[i])
#                 end
#                 nothing
#             end
#         end
#     end
#     for f in (sincos!, tan!, asin!, acos!, atan!, sinhcosh!, tanh!)
#         @eval begin
#             @inline function broadcasted(::$TS{T}, fn::typeof($f), r::$T{T}, a::$T{T}, b::Taylor1{T}, k) where {T}
#                 @inbounds for i in eachindex(k)
#                     fn(r, a, b, k[i])
#                 end
#                 nothing
#             end
#         end
#     end
# end
