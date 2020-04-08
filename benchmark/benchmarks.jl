using BenchmarkTools, TaylorSeries

const SUITE = BenchmarkGroup()

S = SUITE["Constructors"] = BenchmarkGroup()
for i in (05,10,20,50)
    x = rand(i)
    t = tuple(x...)
    S["STaylor1(x::Array[len = $i])"] = @benchmarkable STaylor1($x)
    S["STaylor1(x::Float64, Val($i))"] = @benchmarkable STaylor1($t)
end

function f(g, x, y, n)
    t = x
    for i=1:n
        t = g(x,y)
    end
    t
end
function f(g, x, n)
    t = x
    for i=1:n
        t = g(t)
    end
    t
end
n = 100
dims = (5,20)

S = SUITE["arithmetic"] = BenchmarkGroup()
for i in dims
    r = rand(i)
    x = Taylor1(r)
    x2 = Taylor1(r)
    q = STaylor1(r)
    q2 = STaylor1(r)
    y = rand()
    S["Taylor1{$i,Float64} + Float64"] = @benchmarkable f(+, $x, $y, $n)
    S["STaylor1{$i,Float64} + Float64"] = @benchmarkable f(+, $q, $y, $n)
    S["Taylor1{$i,Float64} - Float64"] = @benchmarkable f(-, $x, $y, $n)
    S["STaylor1{$i,Float64} - Float64"] = @benchmarkable f(-, $q, $y, $n)
    S["Taylor1{$i,Float64} / Float64"] = @benchmarkable f(/, $x, $y, $n)
    S["STaylor1{$i,Float64} / Float64"] = @benchmarkable f(/, $q, $y, $n)
    S["Taylor1{$i,Float64} * Float64"] = @benchmarkable f(*, $x, $y, $n)
    S["STaylor1{$i,Float64} * Float64"] = @benchmarkable f(*, $q, $y, $n)
    S["Taylor1{$i,Float64} + Taylor1{$i,Float64}"] = @benchmarkable f(+, $x, $x2, $n)
    S["STaylor1{$i,Float64} + STaylor1{$i,Float64}"] = @benchmarkable f(+, $q, $q2, $n)
    S["Taylor1{$i,Float64} - Taylor1{$i,Float64}"] = @benchmarkable f(-, $x, $x2, $n)
    S["STaylor1{$i,Float64} - STaylor1{$i,Float64}"] = @benchmarkable f(-, $q, $q2, $n)
    S["Taylor1{$i,Float64} * Taylor1{$i,Float64}"] = @benchmarkable f(*, $x, $x2, $n)
    S["STaylor1{$i,Float64} * STaylor1{$i,Float64}"] = @benchmarkable f(*, $q, $q2, $n)
    S["Taylor1{$i,Float64} / Taylor1{$i,Float64}"] = @benchmarkable f(/, $x, $x2, $n)
    S["STaylor1{$i,Float64} / STaylor1{$i,Float64}"] = @benchmarkable f(/, $q, $q2, $n)
end

S = SUITE["functions"] = BenchmarkGroup()
for i in dims
    r = rand(i)
    x = Taylor1(r)
    q = STaylor1(r)
    y = rand()
    for g in (exp, abs, zero, one, real, imag, conj, adjoint, iszero, isnan,
              isinf, deg2rad, rad2deg)
        S["$(Symbol(g))(Taylor1{Float64}), len = $i"] = @benchmarkable f($g, $x, $n)
        S["$(Symbol(g))(STaylor1{$i,Float64}),"] = @benchmarkable f($g, $q, $n)
    end
end
