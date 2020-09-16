==(a::STaylor1, b::STaylor1) = (a.coeffs == b.coeffs)

iszero(a::STaylor1) = all(iszero, a.coeffs)

zero(::STaylor1{N,T}) where {N, T<:Number} = STaylor1(zero(T), Val(N-1))
one(::STaylor1{N,T}) where {N, T<:Number} = STaylor1(one(T), Val(N-1))

@inline +(a::STaylor1{N,T}, b::STaylor1{N,T}) where {N, T<:Number} = STaylor1(a.coeffs .+ b.coeffs)
@inline -(a::STaylor1{N,T}, b::STaylor1{N,T}) where {N, T<:Number} = STaylor1(a.coeffs .- b.coeffs)
@inline +(a::STaylor1) = a
@inline -(a::STaylor1) = STaylor1(minus_tuple(a.coeffs))

function +(a::STaylor1{N,T}, b::T) where {N, T<:Number}
    STaylor1{N,T}(ntuple(i -> i == 1 ? a.coeffs[1] + b : a.coeffs[i], Val(N)))
end
function +(b::T, a::STaylor1{N,T}) where {N, T<:Number}
    STaylor1{N,T}(ntuple(i -> i == 1 ? a.coeffs[1] + b : a.coeffs[i], Val(N)))
end
function -(a::STaylor1{N,T}, b::T) where {N, T<:Number}
    STaylor1{N,T}(ntuple(i -> i == 1 ? a.coeffs[1] - b : a.coeffs[i], Val(N)))
end
-(b::T, a::STaylor1{N,T}) where {N, T<:Number}  = b + (-a)

#+(a::STaylor1{N,T}, b::S) where {N, T<:NumberNotSeries, S<:NumberNotSeries} = +(promote(a,b)...)
#+(a::STaylor1{N,T}, b::STaylor1{N,S}) where {N, T<:NumberNotSeries, S<:NumberNotSeries} = +(promote(a,b)...)
#+(a::STaylor1{N,T}, b::S) where {N, T<:NumberNotSeries, S<:NumberNotSeries} = +(promote(a,b)...)
#+(b::S, a::STaylor1{N,T}) where {N, T<:NumberNotSeries, S<:NumberNotSeries} = +(promote(b,a)...)

#-(a::STaylor1{N,T}, b::STaylor1{N,S}) where {N, T<:NumberNotSeries, S<:NumberNotSeries} = -(promote(a,b)...)
#-(a::STaylor1{N,T}, b::S) where {N, T<:NumberNotSeries, S<:NumberNotSeries} = -(promote(a,b)...)
#-(b::S, a::STaylor1{N,T}) where {N, T<:NumberNotSeries, S<:NumberNotSeries} = -(promote(b,a)...)

@generated function *(x::STaylor1{N,T}, y::STaylor1{N,T}) where {T<:Number,N}
     ex_calc = quote end
     append!(ex_calc.args, Any[nothing for i in 1:N])
     syms = Symbol[Symbol("c$i") for i in 1:N]
     for j = 1:N
         ex_line = :(x.coeffs[1]*y.coeffs[$j])
         for k = 2:j
             ex_line = :($ex_line + x.coeffs[$k]*y.coeffs[$(j-k+1)])
         end
         sym = syms[j]
         ex_line = :($sym = $ex_line)
         ex_calc.args[j] = ex_line
     end
     exout = :(($(syms[1]),))
     for i = 2:N
         push!(exout.args, syms[i])
     end
     return quote
                Base.@_inline_meta
                $ex_calc
                return STaylor1{N,T}($exout)
             end
end


function *(a::STaylor1{N,T}, b::T) where {N, T<:Number}
    STaylor1{N,T}(b .* a.coeffs)
end
function *(b::T, a::STaylor1{N,T}) where {N, T<:Number}
    STaylor1{N,T}(b .* a.coeffs)
end
