# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries, StaticArrays

using Test
# eeuler = Base.MathConstants.e

@testset "Tests TaylorSeries operations over StaticArrays types" begin
    q = set_variables("q", order=2, numvars=2)
    m = @SMatrix fill(Taylor1(rand(2).*q), 3, 3)
    mt = m'
    @test m isa SMatrix{3, 3, Taylor1{TaylorN{Float64}}, 9}
    @test mt isa SMatrix{3, 3, Taylor1{TaylorN{Float64}}, 9}
    v = @SVector [-1.1, 3.4, 7.62345e-1]
    mtv = mt * v
    @test mtv isa SVector{3, Taylor1{TaylorN{Float64}}}
end
