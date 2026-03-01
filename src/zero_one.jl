# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#


## zero and one ##
for T in (:Taylor1, :HomogeneousPolynomial, :TaylorN)
    @eval iszero(a::$T) = iszero(a.coeffs)
end

function zero(a::Taylor1{T}) where {T<:Number}
    v = Memory{T}(undef, get_order(a)+1)
    v .= zero.(a.coeffs)
    return Taylor1(v)
end
function zero(a::TaylorN{T}) where {T<:NumberNotSeriesN}
    v = Memory{HomogeneousPolynomial{T}}(undef, get_order(a)+1)
    v .= zero.(a.coeffs)
    return TaylorN(v)
end

for T in (:Taylor1, :TaylorN)
    @eval function one(a::$T)
        b = zero(a)
        b[0] = one(b[0])
        return b
    end
end
zero(v::DenseVector{T}) where {T<:AbstractSeries} = zero.(v)

zero(a::HomogeneousPolynomial{T}) where {T<:Number} =
    HomogeneousPolynomial(zero(a.coeffs[1]), get_order(a))

function zeros(a::HomogeneousPolynomial{T}, order::Int) where {T<:Number}
    v = Memory{HomogeneousPolynomial{T}}(undef, order+1)
    z = zero(a[1])
    v .= HomogeneousPolynomial.(z, 0:order)
    return v
end

zeros(::Type{HomogeneousPolynomial{T}}, order::Int) where {T<:Number} =
    zeros( HomogeneousPolynomial(zero(T), 0), order)


one(a::HomogeneousPolynomial{T}) where {T<:Number} =
    HomogeneousPolynomial(one.(a.coeffs), get_order(a))

function ones(a::HomogeneousPolynomial{T}, order::Int) where {T<:Number}
    v = Memory{HomogeneousPolynomial{T}}(undef, order+1)
    for ord in eachindex(v)
        v[ord] = HomogeneousPolynomial(ones(T, size_table[ord]), ord-1)
    end
    return v
end

ones(::Type{HomogeneousPolynomial{T}}, order::Int) where {T<:Number} =
    ones( HomogeneousPolynomial(one(T), 0), order)



# Recursive functions (homogeneous coefficients)
function zero!(a::Taylor1{T}, k::Int) where {T<:NumberNotSeries}
    a[k] = zero(a[k])
    return nothing
end

function zero!(a::Taylor1{T}, k::Int) where {T<:Number}
    for l in eachindex(a[k])
        zero!(a[k], l)
    end
    return nothing
end

function zero!(a::Taylor1{T}) where {T<:Number}
    for k in eachindex(a)
        zero!(a, k)
    end
    return nothing
end


function zero!(a::HomogeneousPolynomial{T}, k::Int) where
        {T<:NumberNotSeries}
    a[k] = zero(a[k])
    return nothing
end
function zero!(a::HomogeneousPolynomial{T}, k::Int) where
        {T<:Number}
    zero!(a[k])
    return nothing
end
function zero!(a::HomogeneousPolynomial{T}) where {T<:Number}
    for k in eachindex(a)
        zero!(a, k)
    end
    return nothing
end


function zero!(a::TaylorN{T}, k::Int) where {T<:Number}
    zero!(a[k])
    return nothing
end

function zero!(a::TaylorN{T}) where {T<:Number}
    for k in eachindex(a)
        zero!(a, k)
    end
    return nothing
end

for T in (:Taylor1, :TaylorN)
    @eval begin
        function zero!(c::$T{T}, a::$T{T}, k::Int) where {T<:Number}
            zero!(c, k)
            return nothing
        end
    end
end


#
function one!(a::Taylor1{T}, k::Int) where {T<:NumberNotSeries}
    a[k] = k == 0 ? one(a[0]) : zero(a[k])
    return nothing
end

function one!(a::Taylor1{T}, k::Int) where {T<:Number}
    zero!(a[k])
    for j in eachindex(a[k])
        one!(a[k], j)
    end
    # one!(a[k])
    return nothing
end

function one!(a::Taylor1{T}) where {T<:Number}
    for k in eachindex(a)
        one!(a, k)
    end
    return nothing
end

function one!(c::TaylorN{T}, k::Int) where {T<:Number}
    zero!(c, k)
    if k == 0
        @inbounds c[0][1] = one(constant_term(c[0][1]))
    end
    return nothing
end


# Taylor1 (including nested Taylor1s)
function one!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where {T<:NumberNotSeries}
    zero!(c, k)
    (k == 0) && (@inbounds c[0] = one(constant_term(a)))
    return nothing
end
function one!(c::Taylor1{T}, a::Taylor1{T}, k::Int) where {T<:Number}
    zero!(c, k)
    k != 0 && return nothing
    for i in eachindex(a[0])
        one!(c[0], a[0], i)
    end
    return nothing
end


function one!(c::TaylorN{T}, a::TaylorN{T}, k::Int) where {T<:Number}
    zero!(c, k)
    (k == 0) && (@inbounds c[0][1] = one(constant_term(a)))
    return nothing
end
