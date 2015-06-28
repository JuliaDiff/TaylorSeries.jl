

facts("Testing an identity proved by Euler (8 variables)") do

    make_variable(name, index::Int) = string(name, TaylorSeries.subscriptify(index))

    variable_names = [make_variable("α", i) for i in 1:4]
    append!(variable_names, [make_variable("β", i) for i in 1:4])

    a1, a2, a3, a4, b1, b2, b3, b4 = set_variables(variable_names, order=4)

    lhs1 = a1^2 + a2^2 + a3^2 + a4^2
    lhs2 = b1^2 + b2^2 + b3^2 + b4^2
    lhs = lhs1 * lhs2

    rhs1 = (a1*b1 - a2*b2 - a3*b3 - a4*b4)^2
    rhs2 = (a1*b2 + a2*b1 + a3*b4 - a4*b3)^2
    rhs3 = (a1*b3 - a2*b4 + a3*b1 + a4*b2)^2
    rhs4 = (a1*b4 + a2*b3 - a3*b2 + a4*b1)^2

    rhs = rhs1 + rhs2 + rhs3 + rhs4

    @fact lhs == rhs  => true
end
