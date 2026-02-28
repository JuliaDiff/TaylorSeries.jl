# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

for T in (:Taylor1, :TaylorN)
    @eval function identity!(c::$T{T}, a::$T{T}) where {T<:Number}
        for k in eachindex(c)
            identity!(c, a, k)
        end
        return nothing
    end
    @eval identity!(c::$T{T}, a::$T{T}, k::Int) where {T<:Number} =
        identity!(c[k], a[k])
end

identity!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where {T<:NumberNotSeries} =
    c[k] = a[k]

@inline function identity!(c::Taylor1{T}, a::T, k::Int) where {T<:NumberNotSeries}
    zero!(c[k])
    k != 0 && return nothing
    c[0] = a
    return nothing
end
function identity!(c::Taylor1{T}, a::T, k::Int) where {T<:Number}
    zero!(c[k])
    k != 0 && return nothing
    for i in eachindex(c[0])
        identity!(c[0], a, i)
    end
    return nothing
end


identity!(c::TaylorN{T}, a::TaylorN{T}, k::Int) where {T<:NumberNotSeries} =
    identity!(c[k], a[k])


function identity!(c::HomogeneousPolynomial{T}, a::HomogeneousPolynomial{T}) where
        {T<:Number}
    for k in eachindex(c)
        identity!(c, a, k)
    end
    return nothing
end

@inline identity!(c::HomogeneousPolynomial{T}, a::HomogeneousPolynomial{T},
    k::Int) where {T<:NumberNotSeries} = c[k] = a[k]

function identity!(c::HomogeneousPolynomial{Taylor1{T}},
        a::HomogeneousPolynomial{Taylor1{T}}, k::Int) where {T<:NumberNotSeries}
    @inbounds for l in eachindex(c[k])
        identity!(c[k], a[k], l)
    end
    return nothing
end
