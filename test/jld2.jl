# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries, JLD2

using Test

function same_taylorn_coefficients(a::TaylorN, b::TaylorN)
    TaylorSeries.order(a) == TaylorSeries.order(b) || return false
    length(a.coeffs) == length(b.coeffs) || return false
    return all(a.coeffs[i].coeffs == b.coeffs[i].coeffs
        for i in eachindex(a.coeffs))
end

@testset "Test TaylorSeries JLD2 extension" begin
    dq = variables!("q", order=4, numvars=6)
    random_TaylorN = [exp(sum(rand(6) .* dq)) for _ in 1:10_000]
    jldsave("test.jld2"; random_TaylorN = random_TaylorN)
    recovered_taylorN = JLD2.load("test.jld2", "random_TaylorN")
    @test recovered_taylorN == random_TaylorN
    rm("test.jld2")
end

@testset "JLD2 explicit JetSpace serialization" begin
    mktempdir() do dir
        path = joinpath(dir, "explicit_spaces.jld2")

        xy_space = JetSpace(order=4, variables=[:x, :y])
        x, y = variables(xy_space)
        f = exp(x + y) + x^2 * y
        g = x - y

        abc_space = JetSpace(order=3, variables=[:a, :b, :c])
        a, b, c = variables(abc_space)
        h = a + b + c

        uv_space_1 = JetSpace(order=2, variables=[:u, :v])
        uv_space_2 = JetSpace(order=2, variables=[:u, :v])
        u1, v1 = variables(uv_space_1)
        u2, v2 = variables(uv_space_2)
        p = u1 + v1
        q = u2 - v2

        jldsave(path; f, g, h, p, q)
        recovered = JLD2.load(path)

        f2 = recovered["f"]
        g2 = recovered["g"]
        h2 = recovered["h"]
        p2 = recovered["p"]
        q2 = recovered["q"]

        @test TaylorSeries.space(f2) === TaylorSeries.space(g2)
        @test TaylorSeries.space(f2) === xy_space
        @test TaylorSeries.get_variable_names(TaylorSeries.space(f2)) == ["x", "y"]
        @test TaylorSeries.order(TaylorSeries.space(f2)) == TaylorSeries.order(xy_space)
        @test same_taylorn_coefficients(f2, f)
        @test same_taylorn_coefficients(g2, g)

        @test TaylorSeries.get_variable_names(TaylorSeries.space(h2)) ==
            ["a", "b", "c"]
        @test TaylorSeries.order(TaylorSeries.space(h2)) ==
            TaylorSeries.order(abc_space)
        @test same_taylorn_coefficients(h2, h)

        @test TaylorSeries.space(p2) !== TaylorSeries.space(q2)
        @test_throws ArgumentError p2 + q2
        @test same_taylorn_coefficients(p2, p)
        @test same_taylorn_coefficients(q2, q)

        ext = Base.get_extension(TaylorSeries, :TaylorSeriesJLD2Ext)
        empty!(getfield(ext, :_write_space_cache))
        empty!(getfield(ext, :_read_space_cache))
        recovered_from_spec = JLD2.load(path)

        f3 = recovered_from_spec["f"]
        g3 = recovered_from_spec["g"]
        h3 = recovered_from_spec["h"]
        p3 = recovered_from_spec["p"]
        q3 = recovered_from_spec["q"]

        @test TaylorSeries.space(f3) === TaylorSeries.space(g3)
        @test TaylorSeries.space(f3) !== xy_space
        @test TaylorSeries.get_variable_names(TaylorSeries.space(h3)) ==
            ["a", "b", "c"]
        @test same_taylorn_coefficients(f3, f)
        @test same_taylorn_coefficients(h3, h)
        @test TaylorSeries.space(p3) !== TaylorSeries.space(q3)
        @test_throws ArgumentError p3 + q3
    end
end
