# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries
using Test
using Random: Xoshiro

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
    tn_res_k = zero(tn_div)
    for k in eachindex(tn_res_k)
        @test isnothing(TS.div!(tn_res_k, tn_div, 2.0, k))
    end
    @test tn_res_k == tn_div / 2.0

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

@testset "Specialized mutating functions" begin
    # The in-place functions and the pair of functions (`s`, `c`) that each one computes
    local SINCOS_INPLACE = ((TaylorSeries.sincos!, sin, cos),
        (TaylorSeries.sincospi!, sinpi, cospi),
        (TaylorSeries.sinhcosh!, sinh, cosh))

    # n-th derivative at 0 of `f`, as a BigFloat
    function deriv_at_zero(f, n::Int)
        f === sin   && return BigFloat((0, 1, 0, -1)[mod(n, 4) + 1])
        f === cos   && return BigFloat((1, 0, -1, 0)[mod(n, 4) + 1])
        f === sinpi && return big(pi)^n * deriv_at_zero(sin, n)
        f === cospi && return big(pi)^n * deriv_at_zero(cos, n)
        f === sinh  && return BigFloat(isodd(n))
        f === cosh  && return BigFloat(iseven(n))
        throw(ArgumentError("unknown function $f"))
    end

    # Coefficient of `t^k * ξ_1^j_1 * ... * ξ_m^j_m` in the expansion of
    # `f(cs * (t + ξ_1 + ... + ξ_m))`: `cs^n f^(n)(0) / (k! j_1! ... j_m!)`, with `n = k + Σ j_i`
    function known_coeff(f, cs, k::Int, js::NTuple{N,Int}) where {N}
        n = k + sum(js; init = 0)
        den = factorial(big(k)) * prod(factorial(big(j)) for j in js; init = big(1))
        return big(cs)^n * deriv_at_zero(f, n) / den
    end

    # All exponent tuples of `nv` variables with total degree ≤ `order`
    monomials(nv::Int, order::Int) =
        [js for js in Iterators.product(ntuple(_ -> 0:order, nv)...) if sum(js) ≤ order]


    @testset "Test sincos!, sincospi!, sinhcosh! for Taylor1{Taylor1}" begin
        # function barrier: `f!` varies in the loop below
        allocs(f!, s, c, a, k) = @allocated f!(s, c, a, k)
        for (f!, fs, fc) in SINCOS_INPLACE, cs in (0.5, 0.5f0, 0.5im),
                (outer, inner) in ((1, 0), (7, 3), (20, 5))
            R = typeof(cs)
            tol = R <: Union{Float32, Complex{Float32}} ? 1e-5 : 1e-13
            u = Taylor1(inner ≥ 1 ? [zero(R), cs] : [zero(R)], inner)
            x = Taylor1([u, Taylor1([cs], inner)], outer)            # cs * (t + u)
            s, c = zero(x), zero(x)
            for k in eachindex(x)
                f!(s, c, x, k)
            end
            @test all(isapprox(s[k][j], known_coeff(fs, cs, k, (j,)); atol = tol)
                for k in 0:outer, j in 0:inner)
            @test all(isapprox(c[k][j], known_coeff(fc, cs, k, (j,)); atol = tol)
                for k in 0:outer, j in 0:inner)
            allocs(f!, s, c, x, outer)                   # warm-up
            @test allocs(f!, s, c, x, outer) == 0
        end
    end

    @testset "Test _mul_scalar2_unchecked!" begin
        sp = JetSpace(order = 4, variables = [:u, :v])
        u, v = variables(sp)
        rng = Xoshiro(20261107)
        randtn() = sum(randn(rng) * u^i * v^j for i in 0:4 for j in 0:4-i)
        function sprinkle_zeros!(x)
            for q in eachindex(x)
                hp = x.coeffs[q+1].coeffs
                r = rand(rng)
                r < 0.25 ? fill!(hp, 0.0) :
                    r < 0.45 ? (hp[rand(rng, eachindex(hp))] = rand(rng, (0.0, -0.0))) : nothing
            end
            return x
        end
        bitequal(x, y) = all(q -> all(map(===, x[q].coeffs, y[q].coeffs)), eachindex(x))

        # `pre * a`, coefficient by coefficient (what the kernel must reproduce with `pre`)
        function prescaled(pre, a)
            pa = zero(a)
            for q in eachindex(a), m in eachindex(a.coeffs[q+1].coeffs)
                pa.coeffs[q+1].coeffs[m] = pre * a.coeffs[q+1].coeffs[m]
            end
            return pa
        end

        for trial in 1:200, scalar in (3, 0, -1.5), pre in (nothing, pi)
            a, b1, b2 = sprinkle_zeros!(randtn()), sprinkle_zeros!(randtn()), sprinkle_zeros!(randtn())
            pa = isnothing(pre) ? a : prescaled(pre, a)
            c1, c2 = randtn(), randtn()
            d1, d2 = zero(c1), zero(c2)              # copies in the same space (no `deepcopy`)
            TS.identity!(d1, c1)
            TS.identity!(d2, c2)
            for k in 0:4
                TS._mul_scalar2_unchecked!(c1, c2, scalar, a, b1, b2, k, pre)
                TS._mul_scalar_unchecked!(d1, scalar, pa, b1, k)
                TS._mul_scalar_unchecked!(d2, scalar, pa, b2, k)
            end
            @test bitequal(c1, d1) && bitequal(c2, d2)
        end
    end

    @testset "Test sincos!, sincospi!, sinhcosh! for Taylor1{TaylorN}" begin
        # function barrier: `f!` varies in the loop below
        allocs(f!, s, c, a, k) = @allocated f!(s, c, a, k)
        for (f!, fs, fc) in SINCOS_INPLACE, cs in (0.5, 0.5im),
                (jetorder, nv) in ((1, 1), (3, 2), (4, 3)), outer in (1, 7, 15)
            sp = JetSpace(order = jetorder, variables = [Symbol(:ξ, i) for i in 1:nv])
            x0 = cs * sum(variables(sp))
            x = Taylor1([x0, cs + zero(x0)], outer)                  # cs * (t + Σ ξ_i)
            x_copy = zero(x)
            for k in eachindex(x)
                TS.identity!(x_copy[k], x[k])
            end
            s, c = zero(x), zero(x)
            for k in eachindex(x)
                f!(s, c, x, k)
            end
            mons = monomials(nv, jetorder)
            @test all(isapprox(getcoeff(s[k], js), known_coeff(fs, cs, k, js); atol = 1e-13)
                for k in 0:outer, js in mons)
            @test all(isapprox(getcoeff(c[k], js), known_coeff(fc, cs, k, js); atol = 1e-13)
                for k in 0:outer, js in mons)
            @test x == x_copy                                        # argument not modified
            allocs(f!, s, c, x, outer)                   # warm-up
            @test allocs(f!, s, c, x, outer) == 0
        end

        # mixed spaces: an error, before anything is modified
        sp1 = JetSpace(order = 2, variables = [:u, :v])
        sp2 = JetSpace(order = 2, variables = [:u, :v])
        u1, _ = variables(sp1)
        u2, _ = variables(sp2)
        for (f!, _, _) in SINCOS_INPLACE
            a1 = Taylor1([0.1 + u1, 0.5 * u1, 0.25 + u1], 2)
            a2 = Taylor1([0.1 + u2, 0.5 * u2, 0.25 + u2], 2)
            s, c = zero(a1), zero(a1)
            f!(s, c, a1, 0)
            @test_throws ArgumentError f!(s, c, a2, 1)
        end
    end

    @testset "sincospi! for TaylorN" begin
        for cs in (0.5, 0.5im), (jetorder, nv) in ((1, 1), (3, 2), (6, 3))
            sp = JetSpace(order = jetorder, variables = [Symbol(:ξ, i) for i in 1:nv])
            x = cs * sum(variables(sp))
            x_copy = zero(x)
            TS.identity!(x_copy, x)
            s1, c1, s2, c2 = zero(x), zero(x), zero(x), zero(x)
            for k in eachindex(x)
                TS.sincospi!(s1, c1, x, k)
            end
            for k in eachindex(x)
                TS.sincospi!(s2, c2, x, k)
            end
            mons = monomials(nv, jetorder)
            @test all(isapprox(getcoeff(s1, js), known_coeff(sinpi, cs, 0, js); atol = 1e-13)
                for js in mons)
            @test all(isapprox(getcoeff(c1, js), known_coeff(cospi, cs, 0, js); atol = 1e-13)
                for js in mons)
            @test x == x_copy                                        # argument not modified
            @test s1 == s2 && c1 == c2                               # same result when repeated
        end
    end
end