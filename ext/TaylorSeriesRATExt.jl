module TaylorSeriesRATExt

using TaylorSeries

isdefined(Base, :get_extension) ? (import RecursiveArrayTools) : (import ..RecursiveArrayTools)

function RecursiveArrayTools.recursivecopy(a::AbstractArray{<:AbstractSeries, N}) where N
    deepcopy(a)
end

end