# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries
using Test

@testset "Mutating functions" begin
    t1 = Taylor1(6)

    # Dictionaries with calls
    @test length(TS._dict_binary_ops) == 5
    # @test length(TS._dict_unary_ops) == 22 # why are these tested?

    @test all([haskey(TS._dict_binary_ops, op)
        for op in [:+, :-, :*, :/, :^]])
    @test all([haskey(TS._dict_binary_calls, op)
        for op in [:+, :-, :*, :/, :^]])

    for kk in keys(TS._dict_binary_ops)
        res = TS._internalmutfunc_call(
            TS._InternalMutFuncs(TS._dict_binary_ops[kk]) )
        @test TS._dict_binary_calls[kk] == res
    end
    for kk in keys(TS._dict_unary_ops)
        res = TS._internalmutfunc_call(
            TS._InternalMutFuncs(TS._dict_unary_ops[kk]) )
        @test TS._dict_unary_calls[kk] == res
    end

    # Some examples
    t1 = Taylor1(5)
    t2 = zero(t1)
    t10 = deepcopy(t1[0])
    TS.pow!(t2, t1, t10, 2, 2)
    @test t2[2] == 1.0
    #
    res = zero(t1)
    TS.add!(res, t1, t2, 3)
    @test res[3] == 0.0
    TS.add!(res, 1, t2, 3)
    @test res[3] == 0.0
    TS.add!(res, t2, 3, 0)
    @test res[0] == 3.0
    TS.subst!(res, t1, t2, 2)
    @test res[2] == -1.0
    TS.subst!(res, t1, 1, 0)
    @test res[0] == -1.0
    @test res[2] == -1.0
    TS.subst!(res, 1, t2, 2)
    @test res[2] == -1.0

    res[3] = rand()
    TS.mul!(res, t1, t2, 3)
    @test res[3] == 1.0
    TS.mul!(res, res, 2, 3)
    @test res[3] == 2.0
    TS.mul!(res, 0.5, res, 3)
    @test res[3] == 1.0

    res[0] = rand()
    TS.div!(res, t2-1, 1+t1, 0)
    res[1] = rand()
    TS.div!(res, t2-1, 1+t1, 1)
    @test res[0] == (t1-1)[0]
    @test res[1] == (t1-1)[1]
    TS.div!(res, res, 2, 0)
    @test res[0] == -0.5

    res = zero(t1)
    TS.identity!(res, t1, 0)
    @test res[0] == t1[0]
    TS.zero!(res, t1, 0)
    TS.zero!(res, t1, 1)
    @test res[0] == zero(t1[0])
    @test res[1] == zero(t1[1])
    TS.one!(res, t1, 0)
    TS.one!(res, t1, 0)
    @test res[0] == one(t1[0])
    @test res[1] == zero(t1[1])

    res = zero(t1)
    TS.abs!(res, -1-t2, 2)
    @test res[2] == 1.0
    @test_throws DomainError TS.abs!(res, t2, 2)

    res = zero(t1)
    TS.abs2!(res, 1-t1, 1)
    @test res[1] == -2.0
    TS.abs2!(res, t1, 2)
    @test res[2] == 1.0

    t2 = Taylor1(Int,15)
    TaylorSeries.zero!(t2)
    @test TaylorSeries.iszero(t2)

end
