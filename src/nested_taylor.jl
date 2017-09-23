# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#


#number of variables in a nested polynomial
function get_nested_numvars{T<:Number}(a::Taylor1{T})
    a_aux = a.coeffs[1]
    dim_a = 1
    while !(typeof(a_aux) <: TaylorSeries.NumberNotSeries)
        dim_a += 1
        a_aux = a_aux.coeffs[1]
    end
    return dim_a
end


function nest{T<:Number}(a::Taylor1{T},num_vars::Integer,dim_a::Integer)
    ord = a.order

    for nest in 1:(num_vars-dim_a)
        a = Taylor1([a],ord)
    end

    return a
end

function nest{T<:Number}(a::Taylor1{T},num_vars::Int64)
    dim_a = get_nested_numvars(a)
    ord = a.order

    for nest in 1:(num_vars-dim_a)
        a = Taylor1([a],ord)
    end

    return a
end

#Shortcut to define an array of independent variables of nested taylors
function set_nested_variables(numvars::Integer,order::Integer,T::Type=Float64)
    nested = Taylor1([zero(T),one(T)],order)
    nested_array = Vector{Number}(numvars)
    numvars_nested = 1

    for i in 1:numvars
        nested_array[i] = nest(nested,numvars,numvars_nested)
        nested = Taylor1([zero(nested),one(nested)],order)
        numvars_nested += 1
    end

    return nested_array
end
