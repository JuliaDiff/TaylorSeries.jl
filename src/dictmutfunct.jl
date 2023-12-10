# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

#=
This file contains some dictionary and functions to build
the dictionaries `_dict_unary_calls` and `_dict_binary_calls`
which allow to call the internal mutating functions. This
may be used to improve memory usage, e.g., to construct the
jet-coefficients used to integrate ODEs.
=#

"""
`_InternalMutFuncs`

Contains parameters and expressions that allow a simple
programmatic construction for calling the internal mutating
functions.
"""
struct _InternalMutFuncs
    namef    :: Symbol   # internal name of the function
    argsf    :: NTuple   # arguments
    defexpr  :: Expr     # defining expr
    auxexpr  :: Expr     # auxiliary expr
end

# Constructor
function _InternalMutFuncs( namef::Tuple )
    if length(namef) == 3
    	return _InternalMutFuncs( namef[1], namef[2], namef[3], Expr(:nothing) )
    else
        return _InternalMutFuncs( namef...)
    end
end


"""
`_dict_binary_ops`

`Dict{Symbol, Array{Any,1}}` with the information to
construct the `_InternalMutFuncs` related to binary
operations.

The keys correspond to the function symbols.

The arguments of the array are the function name (e.g. `add!`), a tuple
with the function arguments, and an `Expr` with the calling pattern. The
convention for the arguments of the functions and the calling pattern
is to use `:_res` for the (mutated) result, `:_arg1` and `_arg2`
for the required arguments, and `:_k` for the computed order
of `:_res`.

"""
const _dict_binary_ops = Dict(
    :+ => (:add!, (:_res, :_arg1, :_arg2, :_k), :(_res = _arg1 + _arg2)),
    :- => (:subst!, (:_res, :_arg1, :_arg2, :_k), :(_res = _arg1 - _arg2)),
    :* => (:mul!, (:_res, :_arg1, :_arg2, :_k), :(_res = _arg1 * _arg2)),
    :/ => (:div!, (:_res, :_arg1, :_arg2, :_k), :(_res = _arg1 / _arg2)),
    :^ => (:pow!, (:_res, :_arg1, :_arg2, :_k), :(_res = _arg1 ^ float(_arg2))),
);

"""
`_dict_binary_ops`

`Dict{Symbol, Array{Any,1}}` with the information to
construct the `_InternalMutFuncs` related to unary
operations.

The keys correspond to the function symbols.

The arguments of the array are the function name (e.g. `add!`), a tuple
with the function arguments, and an `Expr` with the calling pattern. The
convention for the arguments of the functions and the calling pattern
is to use `:_res` for the (mutated) result, `:_arg1`, for the required
argument, possibly `:_aux` when there is an auxiliary expression
needed, and `:_k` for the computed order of `:_res`. When an auxiliary
expression is required, an `Expr` defining its calling pattern is
added as the last entry of the vector.

"""
const _dict_unary_ops = Dict(
    :+ => (:add!,   (:_res, :_arg1, :_k), :(_res = + _arg1)),
    :- => (:subst!, (:_res, :_arg1, :_k), :(_res = - _arg1)),
    :sqr =>  (:sqr!, (:_res, :_arg1, :_k), :(_res = sqr(_arg1))),
    :sqrt => (:sqrt!, (:_res, :_arg1, :_k), :(_res = sqrt(_arg1))),
    :exp =>  (:exp!, (:_res, :_arg1, :_k), :(_res = exp(_arg1))),
    :expm1 =>  (:expm1!, (:_res, :_arg1, :_k), :(_res = expm1(_arg1))),
    :log =>  (:log!, (:_res, :_arg1, :_k), :(_res = log(_arg1))),
    :log1p =>  (:log1p!, (:_res, :_arg1, :_k), :(_res = log1p(_arg1))),
    :identity => (:identity!, (:_res, :_arg1, :_k), :(_res = identity(_arg1))),
    :zero => (:zero!, (:_res, :_arg1, :_k), :(_res = zero(_arg1))),
    :one => (:one!, (:_res, :_arg1, :_k), :(_res = one(_arg1))),
    :abs => (:abs!, (:_res, :_arg1, :_k), :(_res = abs(_arg1))),
    :abs2 => (:abs2!, (:_res, :_arg1, :_k), :(_res = abs2(_arg1))),
    :deg2rad => (:deg2rad!, (:_res, :_arg1, :_k), :(_res = deg2rad(_arg1))),
    :rad2deg => (:rad2deg!, (:_res, :_arg1, :_k), :(_res = rad2deg(_arg1))),
    #
    :sin =>  (:sincos!, (:_res, :_aux, :_arg1, :_k), :(_res = sin(_arg1)),
        :(_aux = cos(_arg1))),
    :cos => (:sincos!, (:_aux, :_res, :_arg1, :_k), :(_res = cos(_arg1)),
        :(_aux = sin(_arg1))),
    :sinpi =>  (:sincospi!, (:_res, :_aux, :_arg1, :_k), :(_res = sinpi(_arg1)),
        :(_aux = cospi(_arg1))),
    :cospi => (:sincospi!, (:_aux, :_res, :_arg1, :_k), :(_res = cospi(_arg1)),
        :(_aux = sinpi(_arg1))),
    :tan => (:tan!, (:_res, :_arg1, :_aux, :_k), :(_res = tan(_arg1)),
        :(_aux = tan(_arg1)^2)),
    :asin => (:asin!, (:_res, :_arg1, :_aux, :_k), :(_res = asin(_arg1)),
        :(_aux = sqrt(1 - _arg1^2))),
    :acos => (:acos!, (:_res, :_arg1, :_aux, :_k), :(_res = acos(_arg1)),
        :(_aux = sqrt(1 - _arg1^2))),
    :atan => (:atan!, (:_res, :_arg1, :_aux, :_k), :(_res = atan(_arg1)),
        :(_aux = 1 + _arg1^2)),
    :sinh => (:sinhcosh!, (:_res, :_aux, :_arg1, :_k), :(_res = sinh(_arg1)),
        :(_aux = cosh(_arg1))),
    :cosh => (:sinhcosh!, (:_aux, :_res, :_arg1, :_k), :(_res = cosh(_arg1)),
        :(_aux = sinh(_arg1))),
    :tanh => (:tanh!, (:_res, :_arg1, :_aux, :_k), :(_res = tanh(_arg1)),
        :(_aux = tanh(_arg1)^2)),
    :asinh => (:asinh!, (:_res, :_arg1, :_aux, :_k), :(_res = asinh(_arg1)),
        :(_aux = sqrt(_arg1^2 + 1))),
    :acosh => (:acosh!, (:_res, :_arg1, :_aux, :_k), :(_res = acosh(_arg1)),
        :(_aux = sqrt(_arg1^2 - 1))),
    :atanh => (:atanh!, (:_res, :_arg1, :_aux, :_k), :(_res = atanh(_arg1)),
        :(_aux = 1 - _arg1^2)),
);



"""
```
_internalmutfunc_call( fn :: _InternalMutFuncs )
```

Creates the appropriate call to the internal mutating
function defined by the `_InternalMutFuncs` object.
This is used to construct [`_dict_unary_calls`](@ref)
and [`_dict_binary_calls`](@ref).
The call contains the prefix `TaylorSeries.`.
"""
_internalmutfunc_call( fn :: _InternalMutFuncs ) = (
    Expr( :call, Meta.parse("TaylorSeries.$(fn.namef)"), fn.argsf... ), fn.defexpr, fn.auxexpr )


"""
`_populate_dicts!()`

Function that populates the internal dictionaries [`_dict_unary_calls`](@ref) and
[`_dict_binary_calls`](@ref)
"""
function _populate_dicts!()
    #Populates the constant vector `_dict_unary_calls`.
    _dict_unary_calls = Dict{Symbol, NTuple{3,Expr}}()
    for kk in keys(_dict_unary_ops)
        res = _internalmutfunc_call( _InternalMutFuncs(_dict_unary_ops[kk]) )
        push!(_dict_unary_calls, kk => res )
    end

    #Populates the constant vector `_dict_binary_calls`.
    _dict_binary_calls = Dict{Symbol, NTuple{3,Expr}}()
    for kk in keys(_dict_binary_ops)
        res = _internalmutfunc_call( _InternalMutFuncs(_dict_binary_ops[kk]) )
        push!(_dict_binary_calls, kk => res )
    end
    return _dict_unary_calls, _dict_binary_calls
end

const _dict_unary_calls, _dict_binary_calls = _populate_dicts!()

@doc """
`_dict_unary_calls::Dict{Symbol, NTuple{2,Expr}}`

Dictionary with the expressions that define the
internal unary functions and the auxiliary functions,
whenever they exist. The keys correspond to those
functions, passed as symbols, with the defined
internal mutating functions.

Evaluating the entries generates expressions that represent
the actual calls to the internal mutating functions.
""" _dict_unary_calls

@doc """
`_dict_binary_calls::Dict{Symbol, NTuple{2,Expr}}`

Dictionary with the expressions that define the
internal binary functions and the auxiliary functions,
whenever they exist. The keys correspond to those
functions, passed as symbols, with the defined
internal mutating functions.

Evaluating the entries generates symbols that represent
the actual calls to the internal mutating functions.
""" _dict_binary_calls

