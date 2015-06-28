# tuple addition:

for n in 1:10
    ex = :( (a[1]+b[1],) )

    for i in 2:n
        push!(ex.args, :(a[$i] + b[$i]))
    end

    f = :( +{T}(a::NTuple{$n,T}, b::NTuple{$n,T}) )

    @eval $f = $ex
end
