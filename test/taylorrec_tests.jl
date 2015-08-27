# testsTaylorRecRec.jl

# using FactCheck
# using TaylorSeries

FactCheck.setstyle(:compact)
# FactCheck.onlystats(true)

x = TaylorRec([0,1])
y = taylorRec_variable(TaylorRec{Int,2},1)
z = @taylorRec(3,2)

facts("\nTaylorRec tests:\nTesting the type and number of variables") do
    @fact get_numvars(y) == 2 --> true
    @fact get_numvars(y.coeffs) == 1 --> true
    @fact get_numvars(x.coeffs) == 0 --> true
    @fact get_numvars(x.coeffs[1]) == 0 --> true
    @fact typeof(x) == TaylorRec{Int,1} --> true
    @fact typeof(y) == TaylorRec{Int,2} --> true
    @fact typeof(z) == TaylorRec{Float64,3} --> true
    @fact eltype(x) == eltype(y) == eltype(y.coeffs) == Int64 --> true
end

facts("Testing TaylorRec constructors and equality") do
    @fact x == @taylorRec(1) --> true
    @fact y == @taylorRec(2,6) --> true
    @fact x == TaylorRec{Int64,1}([0,1]) --> true
    @fact y == TaylorRec{Int64,2}([TaylorRec{Int64,1}([0,1])]) --> true
    @fact x != y --> true
    @fact TaylorRec([1,2,3]) == TaylorRec([1.0,2.0,3.0]) --> true
    @fact TaylorRec{Int,1}(0) == TaylorRec{Int,1}([0]) --> true
    @fact TaylorRec(x) == x --> true
    @fact TaylorRec([x]) == y --> true
    @fact TaylorRec([y]) == TaylorRec(z) --> true
    @fact TaylorRec{Float64,1}(x.coeffs) == taylorRec_variable(TaylorRec{Int,1},1) --> true
    @fact TaylorRec{Float64,3}(z.coeffs) == taylorRec_variable(TaylorRec{Int,3},1) --> true
    @fact_throws AssertionError TaylorRec{Int,0}(0)
    @fact_throws AssertionError TaylorRec{Int,0}([0])
end

facts("Testing zero, one") do
    @fact zero(x) == zero(typeof(x)) == TaylorRec([0,0]) --> true
    @fact one(y) == one(typeof(y)) == TaylorRec([TaylorRec(1)]) --> true
    @fact zero(z) == zero(TaylorRec{Int64,3}) --> true
    @fact one(z) == one(typeof(z)) --> true
end

facts("Testing conversion and promotion") do
    @fact convert(TaylorRec,3) == TaylorRec(3) --> true
    @fact convert(TaylorRec{Float64,1},3) == TaylorRec([3.0]) --> true
    @fact convert(TaylorRec{Float64,2},3) == TaylorRec{Float64,2}(3.0) --> true
    @fact convert(TaylorRec{Float64,2},y) == @taylorRec(2) --> true
    @fact convert(TaylorRec{Int,3},x) == x --> true
    a, b = promote(x,y)
    @fact a == TaylorRec{Int,2}([TaylorRec{Int,1}([0]),TaylorRec{Int,1}([1])]) --> true
    @fact b == TaylorRec{Int,2}([TaylorRec{Int,1}([0,1])]) --> true
    @fact promote(z, 1) == (z, one(z)) --> true
    @fact promote(1, z) == (one(z), z) --> true
end

x = @taylorRec(1, 15, Int)
y = @taylorRec(2, 15, Int)
facts("Testing arithmetics") do
    @fact x == -TaylorRec(-x) == +TaylorRec(x) --> true
    @fact y == -TaylorRec(-y) == +TaylorRec(y) --> true
    @fact x + im == TaylorRec([1im, 1]) --> true
    @fact x - z == TaylorRec{Int64,3}([TaylorRec{Int,2}([TaylorRec{Int,1}([0,-1])]),
        TaylorRec{Int,2}([TaylorRec{Int,1}([1])])]) --> true
    @fact get_numvars(x) == 1 --> true
    @fact 1 + x + x == TaylorRec{Int64,1}([1,2]) --> true
    r = TaylorRec([TaylorRec([TaylorRec([1,1]),TaylorRec([1,0])]),TaylorRec([TaylorRec([1])])])
    @fact 1 + x + y + z == r --> true
    @fact TaylorRec([x+1]) == @taylorRec(2)+1 --> true
    @fact 1-2x == TaylorRec{Int64,1}([1.0,-2.0]) --> true
    @fact get_numvars(x) == 1 --> true
    @fact y*2 == TaylorRec{Int64,2}([TaylorRec([0,2])]) --> true
    @fact get_numvars(y) == 2 --> true
    @fact 2z == TaylorRec([TaylorRec([TaylorRec([0.0,2.0])])]) --> true
    @fact get_numvars(z) == 3 --> true
    xx = @taylorRec(1,4,Int)
    @fact x * xx == TaylorRec([0,0,1]) --> true
    yy = @taylorRec(2,4,Int)
    @fact yy * yy == TaylorRec([TaylorRec([0,0,1])]) --> true
    @fact x * z == TaylorRec{Int64,3}([TaylorRec{Int,2}([TaylorRec{Int,1}([0])]),
        TaylorRec{Int,2}([TaylorRec{Int,1}([0,1])])]) --> true
    @fact z*@taylorRec(3,2) == TaylorRec([TaylorRec([x*xx])]) --> true
    xsquare = x * x
    @fact one(eltype(x)) == x/x --> true
    @fact 1 - x/x == 0 --> true
    @fact (1-xsquare)/(1+x) == 1-x --> true
    yy = y * @taylorRec(2, 15)
    @fact (1-y*y)/(1+y) == 1-y --> true
    @fact 1/(1-y) == TaylorRec([TaylorRec(ones(Int,16))]) --> true
    @fact z*z/z^2 == one(z) --> true
    # @fact_throws AssertionError z/y
    # @fact_throws Error x/y
    @fact (1-x^2)*(1+y^2) ==  1-x^2+y^2-(x*y)^2 --> true
    @fact (x+im)^2 == xsquare + 2*im*x - 1 --> true
    @fact (1im*one(z))^4.0 == 1 --> true
    @fact imag(xsquare+2im*x-1) == 2x --> true
    @fact (x-im)^2 == (xsquare + 2*im*x - 1)' --> true
    @fact real(xsquare + 2*im*x - 1) == xsquare-1 --> true
end

# facts("Testing an identity (proved by Euler)") do
#     for i=1:4
#         ai = symbol(string("a",i))
#         bi = symbol(string("b",i))
#         @eval ($ai) = @taylorRec($i, 4, Int)
#         @eval ($bi) = @taylorRec(4+($i), 4, Int)
#     end
#     expr_lhs1 = a1^2 + a2^2 + a3^2 + a4^2
#     expr_lhs2 = b1^2 + b2^2 + b3^2 + b4^2
#     lhs = expr_lhs1 * expr_lhs2
#     expr_rhs1 = (a1*b1 - a2*b2 - a3*b3 - a4*b4)^2
#     expr_rhs2 = (a1*b2 + a2*b1 + a3*b4 - a4*b3)^2
#     expr_rhs3 = (a1*b3 - a2*b4 + a3*b1 + a4*b2)^2
#     expr_rhs4 = (a1*b4 + a2*b3 - a3*b2 + a4*b1)^2
#     rhs = expr_rhs1 + expr_rhs2 + expr_rhs3 + expr_rhs4
#     @fact lhs == rhs --> true
# end

# exitstatus()
