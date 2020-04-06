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
        t = g(t,x)
        t = g(t,y)
    end
    t
end
function f(g, x, n)
    t = x
    for i=1:n
        t = g(t)
        t = g(t)
    end
    t
end
n = 100

S = SUITE["+"] = BenchmarkGroup()
for i in (5,10,20,50)
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
    S["Taylor1{$i,Float64} + Taylor1{$i,Float64}"] = @benchmarkable f(+, $x, $x2, $n)
    S["STaylor1{$i,Float64} + STaylor1{$i,Float64}"] = @benchmarkable f(+, $q, $q2, $n)
    S["Taylor1{$i,Float64} - Taylor1{$i,Float64}"] = @benchmarkable f(-, $x, $x2, $n)
    S["STaylor1{$i,Float64} - STaylor1{$i,Float64}"] = @benchmarkable f(-, $q, $q2, $n)
end

S = SUITE["functions"] = BenchmarkGroup()
for i in (5,10,20,50)
    r = rand(i)
    x = Taylor1(r)
    q = STaylor1(r)
    y = rand()
    S["Taylor1{$i,Float64}"] = @benchmarkable f(exp, $x, $n)
    S["STaylor1{$i,Float64}"] = @benchmarkable f(exp, $q, $n)
end

#=
S = SUITE["STaylor1 Arithmetic"] = BenchmarkGroup()
for i in (5,10,5)
    S["STaylor1{$i,Float64} Scalar Addition"]
    S["STaylor1{$i,Float64} Scalar Multiplication"]
    S["STaylor1{$i,Float64} Addition"]
    S["STaylor1{$i,Float64} Multiplication"]
end
=#
