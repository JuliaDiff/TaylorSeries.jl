# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Parameters for HomogeneousPolynomial and TaylorN

"Abbreviation for the union of Real and Complex"
const RealOrComplex = Union{Real, Complex}

@doc """
    ParamsTaylorN

DataType holding the current parameters for `TaylorN` and
`HomogeneousPolynomial`.

**Fields:**

- `order          :: Int`  Order (degree) of the polynomials
- `num_vars       :: Int`  Number of variables
- `variable_names :: Array{String,1}` Name of the variables

These parameters can be changed using `set_params_TaylorN(order, numVars)`.
""" ->
type ParamsTaylorN
    order          :: Int
    num_vars       :: Int
    variable_names :: Array{String,1}
end

const _params_TaylorN_ = ParamsTaylorN(6, 2, String["x₁", "x₂"])


## Utilities to get the maximum order and number of variables
get_order() = _params_TaylorN_.order
get_numvars() = _params_TaylorN_.num_vars
get_variable_names() = _params_TaylorN_.variable_names

set_variable_names{T<:AbstractString}(names::Vector{T}) = _params_TaylorN_.variable_names = names

get_variables() = [TaylorN(i) for i in 1:get_numvars()]

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

```jldoctest
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
function set_variables{T<:AbstractString}(R::Type, names::Vector{T}; order=get_order())
    order ≥ 1 || error("Order must be at least 1")

    num_vars = length(names)
    num_vars ≥ 1 || error("Number of variables must be at least 1")

    _params_TaylorN_.variable_names = names

    if !(order == get_order() && num_vars == get_numvars())
        # if these are unchanged, no need to regenerate tables

        _params_TaylorN_.order = order
        _params_TaylorN_.num_vars = num_vars

        resize!(coeff_table,order+1)
        resize!(index_table,order+1)
        resize!(size_table,order+1)
        resize!(pos_table,order+1)

        coeff_table[:], index_table[:], size_table[:], pos_table[:] = generate_tables(num_vars, order)
        gc();
    end

    # return a list of the new variables
    TaylorN{R}[TaylorN(R,i) for i in 1:get_numvars()]
end

set_variables{T}(names::Vector{T}; order=get_order()) =
    set_variables(Float64, names, order=order)

function set_variables{T<:AbstractString}(R::Type, names::T; order=get_order(), numvars=-1)
    variable_names = split(names)

    if length(variable_names) == 1 && numvars ≥ 1
        name = variable_names[1]
        variable_names = [string(name, subscriptify(i)) for i in 1:numvars]
    end

    set_variables(R, variable_names, order=order)
end

set_variables{T<:AbstractString}(names::T; order=get_order(), numvars=-1) =
    set_variables(Float64, names, order=order, numvars=numvars)


@doc """
    show_params_TaylorN()

Display the current parameters for `TaylorN` and `HomogeneousPolynomial` types.
""" ->
function show_params_TaylorN()
    info( """
    Parameters for `TaylorN` and `HomogeneousPolynomial`:
    Maximum order       = $(get_order())
    Number of variables = $(get_numvars())
    Variable names      = $(get_variable_names())
    """)
    nothing
end
