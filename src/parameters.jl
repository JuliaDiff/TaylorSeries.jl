# This file is part of TaylorSeries.jl
#
# Parameters for HomogeneousPolynomial and TaylorN

@doc """Type structure holding the current parameters for `TaylorN` and
`HomogeneousPolynomial`.

    Fieldnames:

    - order:     order (degree) of the polynomials
    - num_vars : number of variables

    These parameters can be changed using `set_params_TaylorN(order,numVars)`
    """ ->

type ParamsTaylorN
    order :: Int
    num_vars  :: Int
    variable_names :: Array{UTF8String,1}
end

const _params_TaylorN_ = ParamsTaylorN(6, 2, UTF8String["x₁", "x₂"])


## Utilities to get the maximum order and number of variables
get_order() = _params_TaylorN_.order
get_numvars() = _params_TaylorN_.num_vars
get_variable_names() = _params_TaylorN_.variable_names

set_variable_names{T<:AbstractString}(names::Vector{T}) = _params_TaylorN_.variable_names = names

get_variables() = [taylorN_variable(i) for i in 1:get_numvars()]

@doc doc"""`set_variables` sets the names and number of the Taylor variables,
as well as the order of the Taylor expansion.""" ->

function set_variables{T<:AbstractString}(R::Type, names::Vector{T}; order=6)
    order >= 1 || error("Order must be at least 1")

    num_vars = length(names)
    num_vars >= 1 || error("Number of variables must be at least 1")

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
    TaylorN{R}[taylorN_variable(R,i) for i in 1:get_numvars()]
end

set_variables{T}(names::Vector{T}; order=6) = set_variables(Float64, names, order=order)

function set_variables{T<:AbstractString}(R::Type, names::T; order=6, numvars=-1)
    variable_names = split(names)

    if length(variable_names) == 1 && numvars >= 1
        name = variable_names[1]
        variable_names = [string(name, subscriptify(i)) for i in 1:numvars]
    end

    set_variables(R, variable_names, order=order)
end

set_variables{T<:AbstractString}(names::T; order=6, numvars=-1) =
    set_variables(Float64, names, order=order, numvars=numvars)

@doc """Display the current parameters for `TaylorN` and
`HomogeneousPolynomial`""" ->
function show_params_TaylorN()
    info( """Parameters for `TaylorN` and `HomogeneousPolynomial`:
    Maximum order       = $(get_order())
    Number of variables = $(get_numvars())
    Variable names      = $(get_variable_names())
    """)
    nothing
end
