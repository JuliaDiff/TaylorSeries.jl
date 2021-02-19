# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries
using Test

@testset "Mutating functions" begin
    t1 = Taylor1(6)

    # Dictionaries with calls
    @test length(TaylorSeries._dict_binary_ops) == 5
    # @test length(TaylorSeries._dict_unary_ops) == 22 # why are these tested?

    @test all([haskey(TaylorSeries._dict_binary_ops, op)
        for op in [:+, :-, :*, :/, :^]])
    @test all([haskey(TaylorSeries._dict_binary_calls, op)
        for op in [:+, :-, :*, :/, :^]])

    for kk in keys(TaylorSeries._dict_binary_ops)
        res = TaylorSeries._internalmutfunc_call(
            TaylorSeries._InternalMutFuncs(TaylorSeries._dict_binary_ops[kk]) )
        @test TaylorSeries._dict_binary_calls[kk] == res
    end
    for kk in keys(TaylorSeries._dict_unary_ops)
        res = TaylorSeries._internalmutfunc_call(
            TaylorSeries._InternalMutFuncs(TaylorSeries._dict_unary_ops[kk]) )
        @test TaylorSeries._dict_unary_calls[kk] == res
    end

    # Some examples
    t1 = Taylor1(5)
    t2 = zero(t1)
    TaylorSeries.pow!(t2, t1, 2, 2)
    @test t2[2] == 1.0
    #
    res = zero(t1)
    TaylorSeries.add!(res, t1, t2, 3)
    @test res[3] == 0.0
    TaylorSeries.add!(res, 1, t2, 3)
    @test res[3] == 0.0
    TaylorSeries.add!(res, t2, 3, 0)
    @test res[0] == 3.0
    TaylorSeries.subst!(res, t1, t2, 2)
    @test res[2] == -1.0
    TaylorSeries.subst!(res, t1, 1, 0)
    @test res[0] == -1.0
    @test res[2] == -1.0
    TaylorSeries.subst!(res, 1, t2, 2)
    @test res[2] == -1.0

    res[3] = rand()
    TaylorSeries.mul!(res, t1, t2, 3)
    @test res[3] == 1.0
    TaylorSeries.mul!(res, res, 2, 3)
    @test res[3] == 2.0
    TaylorSeries.mul!(res, 0.5, res, 3)
    @test res[3] == 1.0

    res[0] = rand()
    TaylorSeries.div!(res, t2-1, 1+t1, 0)
    res[1] = rand()
    TaylorSeries.div!(res, t2-1, 1+t1, 1)
    @test res[0] == (t1-1)[0]
    @test res[1] == (t1-1)[1]
    TaylorSeries.div!(res, res, 2, 0)
    @test res[0] == -0.5

    res = zero(t1)
    TaylorSeries.identity!(res, t1, 0)
    @test res[0] == t1[0]
    TaylorSeries.zero!(res, t1, 0)
    TaylorSeries.zero!(res, t1, 1)
    @test res[0] == zero(t1[0])
    @test res[1] == zero(t1[1])
    TaylorSeries.one!(res, t1, 0)
    TaylorSeries.one!(res, t1, 0)
    @test res[0] == one(t1[0])
    @test res[1] == zero(t1[1])

    res = zero(t1)
    TaylorSeries.abs!(res, -1-t2, 2)
    @test res[2] == 1.0
    @test_throws ArgumentError TaylorSeries.abs!(res, t2, 2)

    res = zero(t1)
    TaylorSeries.abs2!(res, 1-t1, 1)
    @test res[1] == -2.0
    TaylorSeries.abs2!(res, t1, 2)
    @test res[2] == 1.0

end
