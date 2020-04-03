using BenchmarkTools, TaylorSeries

const SUITE = BenchmarkGroup()

S = SUITE["Constructors"] = BenchmarkGroup()
for i in (05,10,50)
    S["STaylor1(x::Array[len = $i])"] = @benchmarkable STaylor1(x) setup = (x = rand($i))
    S["STaylor1(x::Float64, Val($i))"] = @benchmarkable STaylor1(x, Val($i)) setup=(x = rand())
end

function f(g, x, y, n)
    t = x
    for i=1:n
        t = g(t,x)
        t = g(t,y)
    end
    t
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
