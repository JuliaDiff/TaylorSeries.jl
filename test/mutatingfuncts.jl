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
    t1aux = deepcopy(t1)
    TS.pow!(t2, t1, t1aux, 2, 2)
    @test t2[2] == 1.0
    #
    res = zero(t1)
    TS.add!(res, t1, t2)
    @test res == t1 + t2
    TS.subst!(res, t1, t2)
    @test res == t1 - t2
    TS.mul!(res, t1, t2)
    @test res == t1 * t2
    TS.zero!(res)
    TS.muladd!(res, t1, t2)
    @test res == t1 * t2

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

    outer_order = 5
    inner_order = 4
    nt1 = Taylor1([Taylor1([0.5 + i + 0.1*j for j in 0:inner_order],
        inner_order) for i in 0:outer_order], outer_order)
    nt2 = Taylor1([Taylor1([1.0 + 0.2*i - 0.03*j for j in 0:inner_order],
        inner_order) for i in 0:outer_order], outer_order)
    nres = zero(nt1)
    TS.add!(nres, nt1, nt2)
    @test nres == nt1 + nt2
    TS.subst!(nres, nt1, nt2)
    @test nres == nt1 - nt2
    TS.mul!(nres, nt1, nt2)
    @test nres ≈ nt1 * nt2
    TS.zero!(nres)
    for k in eachindex(nres)
        TS.muladd!(nres, nt1, nt2, k)
    end
    @test nres ≈ nt1 * nt2
    ninplace = deepcopy(nt1)
    TS.mul!(ninplace, nt2)
    @test ninplace ≈ nt1 * nt2

    ct1 = 1 + t1 + t1^2
    ct1_tail = collect(ct1[1:order(ct1)])
    @test constant_term!(ct1, 3.5) === ct1
    @test constant_term(ct1) == 3.5
    @test collect(ct1[1:order(ct1)]) == ct1_tail

    nct1 = deepcopy(nt1)
    nct1_tail = [nct1[k] for k in 1:order(nct1)]
    new_inner_constant = Taylor1([2.0, -1.0, 0.25], 2)
    @test constant_term!(nct1, new_inner_constant) === nct1
    @test constant_term(nct1) == new_inner_constant
    @test all(nct1[k] == nct1_tail[k] for k in 1:order(nct1))

    space = JetSpace(order=3, variables=[:x, :y])
    x, y = variables(space)
    tn = 1.0 + x + 2y + 3x*y
    tn_tail = [tn[k] for k in 1:order(tn)]
    @test constant_term!(tn, 7.5) === tn
    @test constant_term(tn) == 7.5
    @test all(tn[k] == tn_tail[k] for k in 1:order(tn))

    hp0 = HomogeneousPolynomial(space, 1.0, 0)
    @test constant_term!(hp0, 4.0) === hp0
    @test hp0[1] == 4.0
    hp1 = HomogeneousPolynomial(space, [1.0, 0.0], 1)
    @test_throws ArgumentError constant_term!(hp1, 2.0)

    tn_div = 2.0 + x + 4y + 6x*y
    tn_res = zero(tn_div)
    @test isnothing(TS.div!(tn_res, tn_div, 2.0))
    @test tn_res == tn_div / 2.0

    t1tn1 = Taylor1([
        (2.0 + 0.1i) + (0.2 + 0.01i)*x + (0.1 - 0.02i)*y +
            (0.03 + 0.01i)*x*y for i in 0:4], 4)
    t1tn2 = Taylor1([
        (1.0 - 0.05i) + (0.3 - 0.02i)*x - (0.04 + 0.01i)*y +
            (0.02 - 0.005i)*x^2 for i in 0:4], 4)
    t1tn_res = zero(t1tn1)
    TS.add!(t1tn_res, t1tn1, t1tn2)
    @test t1tn_res == t1tn1 + t1tn2
    TS.subst!(t1tn_res, t1tn1, t1tn2)
    @test t1tn_res == t1tn1 - t1tn2
    for k in eachindex(t1tn_res)
        TS.mul!(t1tn_res, 2.5, t1tn1, k)
    end
    @test t1tn_res == 2.5 * t1tn1
    for k in eachindex(t1tn_res)
        TS.div!(t1tn_res, t1tn1, 2.5, k)
    end
    t1tn_div_expected = zero(t1tn1)
    for k in eachindex(t1tn_div_expected)
        TS.mul!(t1tn_div_expected, 0.4, t1tn1, k)
    end
    @test t1tn_res ≈ t1tn_div_expected
    d_t1tn = Taylor1(zero(t1tn1[0]), order(t1tn1)-1)
    TS.differentiate!(d_t1tn, t1tn1)
    @test d_t1tn == differentiate(t1tn1)
    log_t1tn = zero(t1tn1)
    for k in eachindex(log_t1tn)
        TS.log!(log_t1tn, t1tn1, k)
    end
    @test exp(log_t1tn) ≈ t1tn1
    small_t1tn = zero(t1tn1)
    for k in eachindex(small_t1tn)
        TS.div!(small_t1tn, t1tn1, 20.0, k)
    end
    log1p_t1tn = zero(small_t1tn)
    for k in eachindex(log1p_t1tn)
        TS.log1p!(log1p_t1tn, small_t1tn, k)
    end
    one_plus_small = zero(small_t1tn)
    for k in eachindex(one_plus_small)
        TS.add!(one_plus_small, 1.0, small_t1tn, k)
    end
    @test exp(log1p_t1tn) ≈ one_plus_small
    @test log1p(small_t1tn) ≈ log1p_t1tn

end
