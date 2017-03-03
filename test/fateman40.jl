# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries
using Base.Test

@testset "Test inspired by Fateman (takes a few seconds)" begin
    x, y, z, w = set_variables(Int128, "x", numvars=4, order=40)

    function fateman2(degree::Int)
        T = Int128
        oneH = HomogeneousPolynomial([one(T)], 0)
        # s = 1 + x + y + z + w
        s = TaylorN(
            [oneH, HomogeneousPolynomial([one(T),one(T),one(T),one(T)],1)], degree)
        s = s^degree
        # s is converted to order 2*ndeg
        s = TaylorN(s, 2*degree)
        return s^2 + s
    end

    function fateman3(degree::Int)
        s = x + y + z + w + 1
        s = s^degree
        s * (s+1)
    end

    f2 = fateman2(20)
    f3 = fateman3(20)
    c = get_coeff(f2,[1,6,7,20])
    @test c == 128358585324486316800
    @test get_coeff(f2,[1,6,7,20]) == c
end
