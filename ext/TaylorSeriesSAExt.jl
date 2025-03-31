module TaylorSeriesSAExt

using TaylorSeries

import Base.promote_op
import LinearAlgebra: matprod

using StaticArrays

promote_op(::typeof(adjoint), ::Type{T}) where {T<:AbstractSeries} = T
promote_op(::typeof(matprod), ::Type{T}, ::Type{U}) where {T <: AbstractSeries, U <: AbstractFloat} = T
promote_op(::typeof(matprod), ::Type{T}, ::Type{T}) where {T <: AbstractSeries} = T

end