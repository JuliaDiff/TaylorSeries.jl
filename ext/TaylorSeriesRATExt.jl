module TaylorSeriesRATExt

using TaylorSeries

import RecursiveArrayTools

function RecursiveArrayTools.recursivecopy(a::AbstractArray{<:AbstractSeries, N}) where N
    deepcopy(a)
end

end