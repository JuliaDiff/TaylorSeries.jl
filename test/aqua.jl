using Test
using TaylorSeries
using Aqua

@testset "Aqua tests (performance)" begin
    # This tests that we don't accidentally run into
    # https://github.com/JuliaLang/julia/issues/29393
    # Aqua.test_unbound_args(TaylorSeries)
    ua = Aqua.detect_unbound_args_recursively(TaylorSeries)
    @test length(ua) == 0

    # See: https://github.com/SciML/OrdinaryDiffEq.jl/issues/1750
    # Test that we're not introducing method ambiguities across deps
    ambs = Aqua.detect_ambiguities(TaylorSeries; recursive = true)
    pkg_match(pkgname, pkdir::Nothing) = false
    pkg_match(pkgname, pkdir::AbstractString) = occursin(pkgname, pkdir)
    filter!(x -> pkg_match("TaylorSeries", pkgdir(last(x).module)), ambs)
    for method_ambiguity in ambs
        @show method_ambiguity
    end
    if VERSION < v"1.10.0-DEV"
        @test length(ambs) == 0
    end
end

@testset "Aqua tests (additional)" begin
    Aqua.test_undefined_exports(TaylorSeries)
    Aqua.test_deps_compat(TaylorSeries)
    Aqua.test_stale_deps(TaylorSeries; ignore=[:Requires])
    Aqua.test_piracies(TaylorSeries)
    Aqua.test_unbound_args(TaylorSeries)
    Aqua.test_project_extras(TaylorSeries)
    Aqua.test_persistent_tasks(TaylorSeries)
end

nothing
