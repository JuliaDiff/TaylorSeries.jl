# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Parameters for HomogeneousPolynomial and TaylorN


"""
    ParamsTaylor1

DataType holding the current variable name(s) for `Taylor1`.
The parameters can be changed using [`set_taylor1_varname`](@ref)
"""
mutable struct ParamsTaylor1
    num_vars :: Int
    var_name :: Vector{String}
end

const _params_Taylor1_ = ParamsTaylor1(1, ["t"])

"""
    set_taylor1_varname(var::String)
    set_taylor1_varname(numvars::Int, names::String)
    set_taylor1_varname(numvars::Int, names::Vector{String})

Change the displayed variable(s) for `Taylor1` objects; the names of
variables cannot include spaces
"""
function set_taylor1_varname(name::String)
    _params_Taylor1_.num_vars = 1
    _params_Taylor1_.var_name = [split(strip(name))[1]]
    return _params_Taylor1_
end

function set_taylor1_varname(numvars::Int, name::String)
    var_names = split(strip(name))
    @assert (length(var_names) == 1) || (length(var_names) == numvars > 0)
    numvars == 1 && return set_taylor1_varname(string(var_names[1]))
    _params_Taylor1_.num_vars = numvars
    if length(var_names) == 1
        _params_Taylor1_.var_name = [var_subscr(string(var_names[1]), i) for i in 1:numvars]
    else
        _params_Taylor1_.var_name = var_names
    end
    return _params_Taylor1_
end

function set_taylor1_varname(numvars::Int, names::Vector{String})
    # _clean_spaces!(names)
    var_names = copy(names)
    for i in eachindex(names)
        var_names[i] = split(strip(var_names[i]))[1]
    end
    @assert length(names) == numvars
    _params_Taylor1_.num_vars = numvars
    _params_Taylor1_.var_name = var_names
    return _params_Taylor1_
end



"""
    ParamsTaylorN

DataType holding the current parameters for `TaylorN` and
`HomogeneousPolynomial`.

**Fields:**

- `order            :: Int`  Order (degree) of the polynomials
- `num_vars         :: Int`  Number of variables
- `variable_names   :: Vector{String}`  Names of the variables
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

"""
    HomogeneousProductTable

Precomputed input-pair and output-major schedules for multiplying two
homogeneous polynomials in a single `TaylorNSpace`.

For multiplying homogeneous polynomials, say, `c = a * b`, each coefficient of
`a` multiplies each coefficient of `b`. So the code walks pairs like::

```
a[1] * b[1]
a[1] * b[2]
a[1] * b[3]
a[2] * b[1]
...
```

For each pair of indices, the field `input_positions[pair]` tells us which output
coefficient `c[pos]` receives that product (in general, many pairs contribute
to the same output coefficient). The input-pair traversing order schedule
maps each pair of input coefficients to the output coefficient it contributes to.
The output-major schedule stores the same products grouped by output coefficient,
so dense kernels (e.g., `_mul_unchecked!`, `_mul_output_major_unchecked!`)
can accumulate one output coefficient at a time.

At initialization, only the input-pair table `input_positions` is filled. Later,
only if a kernel needs the output-major layout, the fields `output_offsets` and
`output_pairs` are filled, i.e., output-major data is filled lazily only when
a kernel requests it. Currently, this lazy loading is realized via a mutable
struct, so that fields can be assigned when appropriate, but it is not
strictly necessary.

# Fields

- `input_positions`: output coefficient position for each input coefficient pair.
- `output_offsets`: start positions of each output coefficient's group in `output_pairs`.
- `output_pairs`: input-pair identifiers grouped by output coefficient.
- `num_right`: number of coefficients in the right input polynomial.
"""
mutable struct HomogeneousProductTable
    input_positions :: Vector{Int}
    output_offsets  :: Vector{Int}
    output_pairs    :: Vector{UInt32}
    num_right       :: Int
end

"""
    TaylorNSpace
    JetSpace

Explicit multivariate Taylor algebra. A space owns the truncation order,
variable metadata, and the lookup tables used by `HomogeneousPolynomial`
and `TaylorN` arithmetic.

# Fields

- `order`: maximum total degree.
- `num_vars`: number of variables.
- `variable_names`: variable names used for printing and compatibility APIs.
- `variable_symbols`: variable symbols corresponding to `variable_names`.
- `coeff_table`: exponent vectors grouped by homogeneous degree.
- `index_table`: hashed exponent-vector labels grouped by homogeneous degree.
- `size_table`: number of monomials for each homogeneous degree.
- `pos_table`: maps hashed exponent labels to coefficient positions.
- `mul_table`: product-table cache indexed by `degree + 1` for the left and right
  homogeneous factors. Entries start as empty placeholders; `input_positions` is
  filled on first product-table use, while `output_offsets` and `output_pairs`
  are filled only if an output-grouped multiplication routine requests them.
- `mul_table_lock`: lock guarding lazy initialization of `mul_table` entries.
"""
mutable struct TaylorNSpace
    order            :: Int
    num_vars         :: Int
    variable_names   :: Vector{String}
    variable_symbols :: Vector{Symbol}
    coeff_table      :: Vector{Vector{Vector{Int}}}
    index_table      :: Vector{Vector{Int}}
    size_table       :: Vector{Int}
    pos_table        :: Vector{Dict{Int,Int}}
    mul_table        :: Vector{Vector{HomogeneousProductTable}}
    mul_table_lock   :: ReentrantLock
end

const JetSpace = TaylorNSpace

const default_space = Ref{TaylorNSpace}()

Base.deepcopy_internal(space::TaylorNSpace, stackdict::IdDict) = space
Base.broadcastable(space::TaylorNSpace) = Ref(space)

function _variable_names_from(names::Vector{T}; numvars::Int=-1) where
        {T<:AbstractString}
    variable_names = String[split(strip(name))[1] for name in names]
    numvars == -1 || length(variable_names) == numvars ||
        error("Number of variable names must match `numvars`")
    return variable_names
end
_variable_names_from(symbs::Vector{T}; numvars::Int=-1) where {T<:Symbol} =
    _variable_names_from(string.(symbs), numvars=numvars)

function _variable_names_from(names::T; numvars::Int=-1) where {T<:AbstractString}
    variable_names = split(names)
    if length(variable_names) == 1 && numvars ≥ 1
        name = variable_names[1]
        variable_names = string.(name, TS.subscriptify.(1:numvars))
    end
    return _variable_names_from(variable_names, numvars=numvars)
end
_variable_names_from(symb::Symbol; numvars::Int=-1) =
    _variable_names_from(string(symb), numvars=numvars)

TaylorNSpace(; order::Int=get_order(), variables) =
    TaylorNSpace(order, _variable_names_from(variables))


## Utilities to get the maximum order, number of variables, their names and symbols
get_order() = _params_TaylorN_.order
get_numvars() = _params_TaylorN_.num_vars
get_variable_names() = _params_TaylorN_.variable_names
get_variable_symbols() = _params_TaylorN_.variable_symbols
get_order(space::TaylorNSpace) = space.order
get_numvars(space::TaylorNSpace) = space.num_vars
get_variable_names(space::TaylorNSpace) = space.variable_names
get_variable_symbols(space::TaylorNSpace) = space.variable_symbols
function lookupvar(s::Symbol)
    ind = findfirst(x -> x==s, _params_TaylorN_.variable_symbols)
    isa(ind, Nothing) && return 0
    return ind
end
function lookupvar(space::TaylorNSpace, s::Symbol)
    ind = findfirst(x -> x==s, space.variable_symbols)
    isa(ind, Nothing) && return 0
    return ind
end


"""
    get_variables([T::Type=Float64], [order::Int=get_order()])
    variables([T::Type=Float64], space::TaylorNSpace; order=get_order(space))
    get_variables([T::Type=Float64], space::TaylorNSpace; order=get_order(space))

Return a `TaylorN{T}` vector with each entry representing an
independent variable. Without an explicit space, `get_variables` uses the
compatibility default space configured by `set_variables`, except that `order`
can be explicitly established by the user without changing internal values for
`num_vars` or `variable_names`.

For explicit spaces, `variables(space)` is the preferred option. The
`get_variables(space)` methods are compatibility aliases for the existing
`get_variables` methods. Omitting `T` defaults to `Float64`.
"""
get_variables(::Type{T}, order::Int=get_order()) where {T} =
    TaylorN.(T, 1:get_numvars(), order=order)
get_variables(order::Int=get_order()) =
    TaylorN.(Float64, 1:get_numvars(), order=order)
variables(::Type{T}, space::TaylorNSpace; order::Int=get_order(space)) where {T} =
    TaylorN.(space, T, 1:get_numvars(space), order=order)
variables(space::TaylorNSpace; order::Int=get_order(space)) =
    variables(Float64, space, order=order)
get_variables(::Type{T}, space::TaylorNSpace; order::Int=get_order(space)) where {T} =
    variables(T, space, order=order)
get_variables(space::TaylorNSpace; order::Int=get_order(space)) =
    variables(space, order=order)

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
3-element Array{TaylorSeries.TaylorN{Int},1}:
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
function set_variables(::Type{R}, names::Vector{T}; order=get_order()) where
        {R, T<:AbstractString}

    space = set_default_space!(TaylorNSpace(order, _variable_names_from(names)))

    # return a list of the new variables
    return variables(R, space)
end
set_variables(::Type{R}, symbs::Vector{T}; order=get_order()) where
    {R,T<:Symbol} = set_variables(R, string.(symbs), order=order)

set_variables(names::Vector{T}; order=get_order()) where {T<:AbstractString} =
    set_variables(Float64, names, order=order)
set_variables(symbs::Vector{T}; order=get_order()) where {T<:Symbol} =
    set_variables(Float64, symbs, order=order)

function set_variables(::Type{R}, names::T; order=get_order(), numvars=-1) where
        {R,T<:AbstractString}

    return set_variables(R, _variable_names_from(names, numvars=numvars), order=order)
end
set_variables(::Type{R}, symbs::Symbol; order=get_order(), numvars=-1) where {R} =
    set_variables(R, string(symbs), order=order, numvars=numvars)

set_variables(names::T; order=get_order(), numvars=-1) where {T<:AbstractString} =
    set_variables(Float64, names, order=order, numvars=numvars)
set_variables(symbs::Symbol; order=get_order(), numvars=-1) =
    set_variables(Float64, string(symbs), order=order, numvars=numvars)


"""
    show_params_TaylorN()

Display the current parameters for `TaylorN` and `HomogeneousPolynomial` types.
"""
show_params_TaylorN() = @info( """
    Parameters for `TaylorN` and `HomogeneousPolynomial`:
    Maximum order       = $(get_order())
    Number of variables = $(get_numvars())
    Variable names      = $(get_variable_names())
    Variable symbols    = $(Symbol.(get_variable_names()))
    """)


# Control the display of the big 𝒪 notation
const bigOnotation = Bool[true]
const _show_default = [false]

"""
    displayBigO(d::Bool) --> nothing

Set/unset displaying of the big 𝒪 notation in  the output
of `Taylor1` and `TaylorN` polynomials. The initial value is
`true`.
"""
displayBigO(d::Bool) = (bigOnotation[end] = d; d)

"""
    use_Base_show(d::Bool) --> nothing

Use `Base.show_default` method (default `show` method
in Base), or a custom display. The initial value is
`false`, so customized display is used.
"""
use_show_default(d::Bool) = (_show_default[end] = d; d)
