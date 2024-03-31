# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries, JLD2

using Test

@testset "Test TaylorSeries JLD2 extension" begin
    dq = set_variables("q", order=4, numvars=6)
    random_TaylorN = [cos(sum(dq .* rand(6))), sin(sum(dq .* rand(6))), tan(sum(dq .* rand(6)))]
    jldsave("test.jld2"; random_TaylorN = random_TaylorN)
    recovered_taylorN = JLD2.load("test.jld2", "random_TaylorN")
    @test recovered_taylorN == random_TaylorN
    rm("test.jld2")
end
