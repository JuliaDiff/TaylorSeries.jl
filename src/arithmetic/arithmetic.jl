# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

include("add_subst.jl")
include("multiplication.jl")
include("division.jl")
include("sqr.jl")
include("sqrt.jl")
include("power.jl")


"""
    mul!(Y, A, B)

Multiply A*B and save the result in Y.
"""
function mul!(y::Vector{Taylor1{T}},
        a::Union{Matrix{T},SparseMatrixCSC{T}},
        b::Vector{Taylor1{T}}) where {T<:Number}

    n, k = size(a)
    @assert (length(y)== n && length(b)== k)

    # determine the maximal order of b
    order = maximum(TS.order.(b))

    # Use matrices of coefficients (of proper size) and mul!
    # B = zeros(T, k, order+1)
    B = Array{T}(undef, k, order+1)
    B = zero.(B)
    for i = 1:k
        @inbounds ord = TS.order(b[i])
        @inbounds for j = 1:ord+1
            B[i,j] = b[i][j-1]
        end
    end
    Y = Array{T}(undef, n, order+1)
    mul!(Y, a, B)
    @inbounds for i = 1:n
        # y[i] = Taylor1( collect(Y[i,:]), order)
        y[i] = Taylor1( Y[i,:], order)
    end

    return y
end


# Adapted from (Julia v1.2) stdlib/v1.2/LinearAlgebra/src/dense.jl#721-734,
# licensed under MIT "Expat".
# Specialize a method of `inv` for Matrix{Taylor1{T}}. Simply, avoid pivoting,
# since the polynomial field is not an ordered one.
# function Base.inv(A::StridedMatrix{Taylor1{T}}) where T
#     checksquare(A)
#     S = Taylor1{typeof((one(T)*zero(T) + one(T)*zero(T))/one(T))}
#     AA = convert(AbstractArray{S}, A)
#     if istriu(AA)
#         Ai = triu!(parent(inv(UpperTriangular(AA))))
#     elseif istril(AA)
#         Ai = tril!(parent(inv(LowerTriangular(AA))))
#     else
#         # Do not use pivoting !!
#         Ai = inv!(lu(AA, Val(false)))
#         Ai = convert(typeof(parent(Ai)), Ai)
#     end
#     return Ai
# end

# see https://github.com/JuliaLang/julia/pull/40623
const LU_RowMaximum = RowMaximum()
const LU_NoPivot = NoPivot()

# Adapted from (Julia v1.2) stdlib/v1.2/LinearAlgebra/src/lu.jl#240-253
# and (Julia v1.4.0-dev) stdlib/LinearAlgebra/v1.4/src/lu.jl#270-274,
# licensed under MIT "Expat".
# Specialize a method of `lu` for Matrix{Taylor1{T}}, which avoids pivoting,
# since the polynomial field is not an ordered one.
# We can't assume an ordered field so we first try without pivoting
function lu(A::AbstractMatrix{Taylor1{T}}; check::Bool = true) where {T<:Number}
    S = Taylor1{lutype(T)}
    F = lu!(copy_oftype(A, S), LU_NoPivot; check = false)
    if issuccess(F)
        return F
    else
        return lu!(copy_oftype(A, S), LU_RowMaximum; check = check)
    end
end


# Fast allocation-free matrix multiplication
for T in (:Taylor1, :TaylorN)
    @eval function matmul!(C::Matrix{$T{T}},
                           A::Matrix{$T{T}}, B::Matrix{$T{T}}) where {T}
        mc, nc = size(C)
        ma, na = size(A)
        mb, nb = size(B)
        @assert (na == mb && mc == ma && nc == nb)
        for j in axes(C,2)
            for i in axes(C,1)
                TS.zero!(C[i,j])
                for k in 1:na
                    TS.muladd!(C[i,j], A[i,k], B[k,j])
                end
            end
        end
        return nothing
    end
end
