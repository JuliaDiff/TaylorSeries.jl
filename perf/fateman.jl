using TaylorSeries

set_variables("x", numvars=4, order=40)

function fateman1(degree::Int)
    T = Int128
    oneH = HomogeneousPolynomial(one(T), 0)
    # s = 1 + x + y + z + w
    s = TaylorN( [oneH, HomogeneousPolynomial([one(T),one(T),one(T),one(T)],1)], degree )
    s = s^degree
    # s is converted to order 2*ndeg
    s = TaylorN(s, 2*degree)

    s * ( s+TaylorN(oneH, 2*degree) )
end

f1 = fateman1(0)
println("Fateman 1:")
@time f1 = fateman1(20)

function fateman2(degree::Int)
    T = Int128
    oneH = HomogeneousPolynomial(one(T), 0)
    # s = 1 + x + y + z + w
    s = TaylorN( [oneH, HomogeneousPolynomial([one(T),one(T),one(T),one(T)],1)], degree )
    s = s^degree
    # s is converted to order 2*ndeg
    s = TaylorN(s, 2*degree)
    return s^2 + s
end

f2 = fateman2(0);
println("Fateman 2:")
@time f2 = fateman2(20);
