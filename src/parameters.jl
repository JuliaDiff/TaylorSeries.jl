# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Parameters for HomogeneousPolynomial and TaylorN


"""
    ParamsTaylorN

DataType holding the current parameters for `TaylorN` and
`HomogeneousPolynomial`.

**Fields:**

- `order            :: Int`  Order (degree) of the polynomials
- `num_vars         :: Int`  Number of variables
- `variable_names   :: Vector{String}` Names of the variables
- `variable_symbols :: Vector{Symbol}`  Symbols of the variables

These parameters can be changed using [`set_variables`](@ref)
"""
mutable struct ParamsTaylorN
    order            :: Int
    num_vars         :: Int
    variable_names   :: Vector{String}
    variable_symbols :: Vector{Symbol}
end


ParamsTaylorN(order, num_vars, variable_names) = ParamsTaylorN(order, num_vars, variable_names, Symbol.(variable_names))

const _params_TaylorN_ = ParamsTaylorN(6, 2, ["x₁", "x₂"])


## Utilities to get the maximum order, number of variables, their names and symbols
get_order() = _params_TaylorN_.order
get_numvars() = _params_TaylorN_.num_vars
get_variable_names() = _params_TaylorN_.variable_names
get_variable_symbols() = _params_TaylorN_.variable_symbols
function lookupvar(s::Symbol)
    @compat ind = findfirst(x -> x==s, _params_TaylorN_.variable_symbols)
    @compat isa(ind, Nothing) && return 0
    return ind
end

function set_variable_names(varnames::Vector{T}) where {T<:AbstractString}
    _params_TaylorN_.variable_names = varnames
    _params_TaylorN_.variable_symbols = Symbol.(varnames)
    nothing
end
"""
    get_variables(T::Type, [order::Int=get_order()])

Return a `TaylorN{T}` vector with each entry representing an
independent variable. It takes the default `_params_TaylorN_` values
if `set_variables` hasn't been changed with the exception that `order`
can be explicitely established by the user without changing internal values
for `num_vars` or `variable_names`. Ommiting `T` defaults to `Float64`.
"""
get_variables(T::Type, order::Int=get_order()) =
    [TaylorN(T, i, order=order) for i in 1:get_numvars()]
get_variables(order::Int=get_order()) =
    [TaylorN(Float64, i, order=order) for i in 1:get_numvars()]

"""
    set_variables([T::Type], names::String; [order=get_order(), numvars=-1])

Return a `TaylorN{T}` vector with each entry representing an
independent variable. `names` defines the output for each variable
(separated by a space). The default type `T` is `Float64`,
and the default for `order` is the one defined globally.
Changing the `order` or `numvars` resets the hash_tables.

If `numvars` is not specified, it is inferred from `names`. If only
one variable name is defined and `numvars>1`, it uses this name with
subscripts for the different variables.

```julia
julia> set_variables(Int, "x y z", order=4)
3-element Array{TaylorSeries.TaylorN{Int64},1}:
  1 x + 𝒪(‖x‖⁵)
  1 y + 𝒪(‖x‖⁵)
  1 z + 𝒪(‖x‖⁵)

julia> set_variables("α", numvars=2)
2-element Array{TaylorSeries.TaylorN{Float64},1}:
  1.0 α₁ + 𝒪(‖x‖⁵)
  1.0 α₂ + 𝒪(‖x‖⁵)

julia> set_variables("x", order=6, numvars=2)
2-element Array{TaylorSeries.TaylorN{Float64},1}:
  1.0 x₁ + 𝒪(‖x‖⁷)
  1.0 x₂ + 𝒪(‖x‖⁷)
```
"""
function set_variables(R::Type, names::Vector{T}; order=get_order()) where
        {T<:AbstractString}

    order ≥ 1 || error("Order must be at least 1")

    num_vars = length(names)
    num_vars ≥ 1 || error("Number of variables must be at least 1")

    _params_TaylorN_.variable_names = names
    _params_TaylorN_.variable_symbols = Symbol.(names)


    if !(order == get_order() && num_vars == get_numvars())
        # if these are unchanged, no need to regenerate tables

        _params_TaylorN_.order = order
        _params_TaylorN_.num_vars = num_vars

        resize!(coeff_table,order+1)
        resize!(index_table,order+1)
        resize!(size_table,order+1)
        resize!(pos_table,order+1)

        coeff_table[:], index_table[:], size_table[:], pos_table[:] =
            generate_tables(num_vars, order)
        @compat GC.gc();
    end

    # return a list of the new variables
    TaylorN{R}[TaylorN(R,i) for i in 1:get_numvars()]
end
set_variables(R::Type, symbs::Vector{T}; order=get_order()) where
    {T<:Symbol} = set_variables(R, string.(symbs), order=order)

set_variables(names::Vector{T}; order=get_order()) where {T<:AbstractString} =
    set_variables(Float64, names, order=order)
set_variables(symbs::Vector{T}; order=get_order()) where {T<:Symbol} =
    set_variables(Float64, symbs, order=order)

function set_variables(R::Type, names::T; order=get_order(), numvars=-1) where
        {T<:AbstractString}

    variable_names = split(names)

    if length(variable_names) == 1 && numvars ≥ 1
        name = variable_names[1]
        variable_names = [string(name, subscriptify(i)) for i in 1:numvars]
    end

    set_variables(R, variable_names, order=order)
end
set_variables(R::Type, symbs::Symbol; order=get_order(), numvars=-1) =
    set_variables(R, string(symbs), order=order, numvars=numvars)

set_variables(names::T; order=get_order(), numvars=-1) where {T<:AbstractString} =
    set_variables(Float64, names, order=order, numvars=numvars)
set_variables(symbs::Symbol; order=get_order(), numvars=-1) =
    set_variables(Float64, string(symbs), order=order, numvars=numvars)


"""
    show_params_TaylorN()

Display the current parameters for `TaylorN` and `HomogeneousPolynomial` types.
"""
function show_params_TaylorN()
    Compat.@info( """
    Parameters for `TaylorN` and `HomogeneousPolynomial`:
    Maximum order       = $(get_order())
    Number of variables = $(get_numvars())
    Variable names      = $(get_variable_names())
    Variable symbols    = $(Symbol.(get_variable_names()))
    """)
    nothing
end


# Control the display of the big 𝒪 notation
const bigOnotation = Bool[true]

"""
    displayBigO(d::Bool) --> nothing

Set/unset displaying of the big 𝒪 notation in  the output
of `Taylor1` and `TaylorN` polynomials. The initial value is
`true`.
"""
displayBigO(d::Bool) = (bigOnotation[end] = d; d)
