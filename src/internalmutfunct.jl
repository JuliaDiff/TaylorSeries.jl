"""
`_InternalMutFuncs`

Contains parameters and expressions that allow a simple
programatic construction for calling the internal mutating
functions.
"""
struct _InternalMutFuncs
    namef    :: Symbol     # internal name of the function
    argsf    :: NTuple   # arguments
    auxbool  :: Bool     # Needs auxiliary Taylor object?
	auxpos   :: Int      # Position in the arg list of the aux Taylor object
    auxdef   :: Expr     # Definition of the aux Taylor object
    # compbool :: Bool     # Component of the result?
    # indx     :: Int      # Actual component
end

# Constructors
function _InternalMutFuncs( namef::Vector )
    if length(namef) == 3
    	return _InternalMutFuncs( namef[1], namef[2], false, 0, Expr(:nothing) )
    else
        return _InternalMutFuncs( namef[1], namef[2], namef[end]...)
    end
end


"""
`__dict_internalmutfuncs`

`Dict{Symbol, Array{Any,1}}` with the information to
construct the different `_InternalMutFuncs` related to binary
operations.

The keys correspond to the function symbols.

The arguments of the array are the internal call (function
name and arguments) for the internal mutating function defined,
an expression with the calling pattern related of the function
and, in case when the corresponding `_InternalMutFuncs` object
has either `auxbool` or `compbool` `true`, the remaining
arguments to define the `_InternalMutFuncs` object.
"""
const _dict_internalmutfuncs = Dict(
    :+ => [:add!, (:_res, :_arg1, :_arg2, :_k), :(_res = _arg1 + _arg2)],
    :- => [:subst!, (:_res, :_arg1, :_arg2, :_k), :(_res = _arg1 - _arg2)],
    :* => [:mul!, (:_res, :_arg1, :_arg2, :_k), :(_res = _arg1 * _arg2)],
    :/ => [:div!, (:_res, :_arg1, :_arg2, :_k), :(_res = _arg1 / _arg2)],
    :^ => [:pow!, (:_res, :_arg1, :_arg2, :_k), :(_res = _arg1 ^ _arg2)],
    :sqr =>  [:sqr!, (:_res, :_arg1, :_k), :(_res = sqr(_arg1))],
    :sqrt => [:sqrt!, (:_res, :_arg1, :_k), :(_res = sqrt(_arg1))],
    :exp =>  [:exp!, (:_res, :_arg1, :_k), :(_res = exp(_arg1))],
    :log =>  [:log!, (:_res, :_arg1, :_k), :(_res = log(_arg1))],
    :sin =>  [:sincos!, (:_res, :_aux, :_arg1, :_k), :(_res = sin(_arg1)),
        (true, 2, :(_aux = cos(constant_term(:_arg1))))],
    # :cos => [:sincos!, (:_aux, :_res, :_arg1, :_k), :(_res = cos(_arg1)),
    #     (true, 1, :(_aux = sin(constant_term(_arg1))), true, 2)],
    # :tan => [:tan!, (:_res, :_arg1, :_aux, :_k), :(_res = tan(_arg1)),
    #     (true, 3, :(_aux = tan(constant_term(_arg1))^2), false, 1)],
    # :asin => [:asin!, (:_res, :_arg1, :_aux, :_k), :(_res = asin(_arg1)),
    #     (true, 3, :(_aux = sqrt(1 - constant_term(_arg1)^2)), false, 1)],
    # :acos => [:acos!, (:_res, :_arg1, :_aux, :_k), :(_res = acos(_arg1)),
    #     (true, 3, :(_aux = sqrt(1 - constant_term(_arg1)^2)), false, 1)],
    # :atan => [:atan!, (:_res, :_arg1, :_aux, :_k), :(_res = atan(_arg1)),
    #     (true, 3, :(_aux = 1 + constant_term(_arg1)^2), false, 1)],
    # :sinh => [:sincosh!, (:_res, :_aux, :_arg1, :_k), :(_res = sinh(_arg1)),
    #     (true, 2, :(_aux = cosh(constant_term(_arg1))), true, 1)],
    # :cosh => [:sincosh!, (:_aux, :_res, :_arg1, :_k), :(_res = cosh(_arg1)),
    #     (true, 1, :(_aux = sinh(constant_term(_arg1))), true, 2)],
    # :tanh => [:tanh!, (:_res, :_arg1, :_aux, :_k), :(_res = tanh(_arg1)),
    #     (true, 3, :(_aux = tanh( constant_term(_arg1) )^2), true, 2)],
);



"""
```
_internalmutfunc_call( fn :: _InternalMutFuncs )
```

Creates the appropriate call to the internal mutating
function defined by the `_InternalMutFuncs` object.
This is used to construct
[`_dict_internalcalls`](@ref).
"""
function _internalmutfunc_call( fn :: _InternalMutFuncs )
    # ex = ifelse(~fn.compbool,
    #     Expr( :call, fn.namef, fn.argsf... ),
	# 	Expr( :ref, Expr(:call, fn.namef, fn.argsf... ), fn.indx ) )

    ex = Expr( :call, fn.namef, fn.argsf... )
	fn.auxbool && @assert(fn.auxdef.args[1] == fn.argsf[fn.auxpos])
	exaux = fn.auxdef

	return ex, exaux
end



"""
`_dict_internalcalls::Dict{Symbol, NTuple{2,Expr}}`

Dictionary with the expressions that define the
internal functions and the auxiliary functions, whenever
they exist. The keys correspond to those functions, passed
as symbols, with defined internal mutating functions.

Evaluating
the entries generates symbols that represent the
actual call to the internal mutating functions.
"""
_dict_internalcalls = Dict{Symbol, NTuple{2,Expr}}()

#Populates the constant vector `_dict_internalcalls`.
for kk in keys(_dict_internalmutfuncs)
    res = _internalmutfunc_call( _InternalMutFuncs(_dict_internalmutfuncs[kk]) )
    push!(_dict_internalcalls, kk => res )
end
