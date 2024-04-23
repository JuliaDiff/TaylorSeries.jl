# This file is part of TaylorSeries.jl, MIT licensed
#

using TaylorSeries, RecursiveArrayTools

using Test

@testset "Tests TaylorSeries RecursiveArrayTools extension" begin
    dq = get_variables()
    x = Taylor1([0.9+2dq[1],-1.1dq[1], 0.7dq[2], 0.5dq[1]-0.45dq[2],0.9dq[1]])
    xx = [x,x]
    yy = recursivecopy(xx)
    @test yy == xx # yy and xx are equal...
    @test yy !== xx # ...but they're not the same object in memory
end
