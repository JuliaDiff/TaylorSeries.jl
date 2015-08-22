using TaylorSeries

const order = 20
const x, y, z, w = set_variables(Int128, "x", numvars=4, order=2order)

function fateman1(degree::Int)
    T = Int128
    oneH = HomogeneousPolynomial(one(T), 0)
    # s = 1 + x + y + z + w
    s = TaylorN( [oneH, HomogeneousPolynomial([one(T),one(T),one(T),one(T)],1)], degree )
    s = s^degree
    # s is converted to order 2*ndeg
    s = TaylorN(s, 2*degree)

    s * ( s+TaylorN(oneH, 2*degree) )
end

function fateman2(degree::Int)
    T = Int128
    oneH = HomogeneousPolynomial(one(T), 0)
    # s = 1 + x + y + z + w
    s = TaylorN( [oneH, HomogeneousPolynomial([one(T),one(T),one(T),one(T)],1)], degree )
    s = s^degree
    # s is converted to order 2*ndeg
    s = TaylorN(s, 2*degree)
    return s^2 + s
end

function fateman3(degree::Int)
    s = x + y + z + w + 1
    s = s^degree
    s * (s+1)
end

function fateman4(degree::Int)
    s = x + y + z + w + 1
    s = s^degree
    s^2 + s
end

function run_fateman(N)
    results = Any[]
    nn = 5
    for f in (fateman1, fateman2, fateman3, fateman4)
        f(0)
        println("Running $f")
        @time result = f(N)
        # push!(results, result) # This may take a lot of memory
        t = Inf
        tav = 0.0
        for i = 1:nn
            ti = @elapsed f(N)
            tav += ti
            t = min(t,ti)
        end
        println("\tAverage time of $nn runs: ", tav/nn)
        println("\tMinimum time of $nn runs: ", t)
    end
    results
end


println("Running Fateman with order $order...")

results = run_fateman(order);

println("Done.")

# @assert results[1] == results[2] == results[3] == results[4]
