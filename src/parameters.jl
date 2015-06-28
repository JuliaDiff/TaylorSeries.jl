@doc """Type structure holding the current parameters for `TaylorN` and
`HomogeneousPolynomial`.

    Fieldnames:

    - order: maximum order (degree) of the polynomials
    - num_vars : maximum number of variables

    These parameters can be changed using `set_params_TaylorN(order,num_vars)`
    """ ->
type ParamsTaylorN
    order :: Int
    num_vars  :: Int
    variable_names :: Array{UTF8String,1}

    ParamsTaylorN() = new()  # empty constructor
end

@doc """Display the current parameters for `TaylorN` and
`HomogeneousPolynomial`""" ->
function show_params_TaylorN()
    info( """Parameters for `TaylorN` and `HomogeneousPolynomial`:
    Maximum order       = $(_params_TaylorN_.order)
    Number of variables = $(_params_TaylorN_.num_vars)
    Variable names      = $(_params_TaylorN_.variable_names)
    """)
    nothing
end

global const _params_TaylorN_ = ParamsTaylorN()  # 6, 2, UTF8String["x₁", "x₂"])


get_numvars() = _params_TaylorN_.num_vars
get_order() = _params_TaylorN_.order

get_variable_names() = _params_TaylorN_.variable_names
set_variable_names{T<:String}(names::Vector{T}) = _params_TaylorN_.variable_names = names


@doc doc"""`set_variables` sets the names and number of the Taylor variables,
as well as the order of the Taylor expansion.""" ->

function set_variables{T}(names::Vector{T}; order=6)

    global index_table, size_table, pos_table

    order >= 1 || error("Order must be at least 1")

    num_vars = length(names)
    num_vars >= 1 || error("Number of variables must be at least 1")

    _params_TaylorN_.variable_names = names

    if !(order == _params_TaylorN_.order && num_vars == _params_TaylorN_.num_vars)
        # if these are unchanged, no need to regenerate tables

        _params_TaylorN_.order = order
        _params_TaylorN_.num_vars = num_vars


        index_table, size_table, pos_table = generate_tables(num_vars, order)
        gc();
    end

    # return a list of the new variables
    [taylorN_variable(i) for i in 1:get_numvars()]
end


function set_variables{T<:String}(names::T; order=6, numvars=-1)
    variable_names = split(names)

    if length(variable_names) == 1 && numvars >= 1
        name = variable_names[1]
        variable_names = [string(name, subscriptify(i)) for i in 1:numvars]
    end

    set_variables(variable_names, order=order)

end

