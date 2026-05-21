# Product benchmarks for TaylorN and nested Taylor types.
#
# Run from the repository root with:
#     julia --project=. perf/product_benchmarks.jl
#
# This script requires BenchmarkTools to be available in the active or shared
# Julia environment. It intentionally does not add BenchmarkTools as a package
# dependency of TaylorSeries.jl.

try
    import BenchmarkTools
catch err
    if err isa ArgumentError
        error("""
        BenchmarkTools is required to run this benchmark.

        Install it in a shared or temporary environment, for example:
            julia -e 'using Pkg; Pkg.add("BenchmarkTools")'

        Then rerun:
            julia --project=. perf/product_benchmarks.jl
        """)
    end
    rethrow()
end

using TaylorSeries

const DA_ORDER = parse(Int, get(ENV, "TS_PRODUCT_BENCH_DA_ORDER", "6"))
const TIME_ORDER = parse(Int, get(ENV, "TS_PRODUCT_BENCH_TIME_ORDER", "7"))
const SECONDS = parse(Float64, get(ENV, "TS_PRODUCT_BENCH_SECONDS", "1.0"))
const SAMPLES = parse(Int, get(ENV, "TS_PRODUCT_BENCH_SAMPLES", "10000"))

function lift_taylor1(a::TaylorN{Float64}, orderT::Int)
    coeffs = Vector{HomogeneousPolynomial{Taylor1{Float64}}}(undef, get_order(a)+1)
    for ord in eachindex(a)
        hp = a[ord]
        hpc = Vector{Taylor1{Float64}}(undef, length(hp))
        for i in eachindex(hp)
            c = hp[i]
            hpc[i] = Taylor1([c, (0.05 + 0.01*ord) * c], orderT)
        end
        coeffs[ord+1] = HomogeneousPolynomial(a.space, hpc, ord)
    end
    return TaylorN(a.space, coeffs, get_order(a))
end

function product_inputs(; da_order::Int=DA_ORDER, time_order::Int=TIME_ORDER)
    space = JetSpace(order=da_order, variables=[:x1, :x2, :x3, :x4])
    x = variables(space)

    f = sin(x[1] + x[2]*x[3]) + exp(x[4]) + (x[1]-x[2]+x[3])^3
    g = cos(x[2] - x[4]) + x[1]^2 - x[3]*x[4] + exp(x[1]*x[2])
    resN = zero(f)

    tf = Taylor1([f + (0.03*i) * g for i in 0:time_order], time_order)
    tg = Taylor1([g - (0.02*i) * f for i in 0:time_order], time_order)
    resTn = zero(tf)

    fn_t1 = lift_taylor1(f, time_order)
    gn_t1 = lift_taylor1(g, time_order)
    resNt1 = zero(fn_t1)

    return (; f, g, resN, tf, tg, resTn, fn_t1, gn_t1, resNt1)
end

function report(label::AbstractString, trial)
    min_est = minimum(trial)
    med_est = BenchmarkTools.median(trial)
    println(rpad(label, 42),
        " min=", BenchmarkTools.prettytime(min_est.time),
        " median=", BenchmarkTools.prettytime(med_est.time),
        " memory=", BenchmarkTools.prettymemory(med_est.memory),
        " allocs=", med_est.allocs)
    return nothing
end

function bench_product(label::AbstractString, a, b)
    trial = BenchmarkTools.@benchmark $a * $b evals=1 seconds=SECONDS samples=SAMPLES
    report(label, trial)
    return trial
end

function bench_mul!(label::AbstractString, res, a, b)
    trial = BenchmarkTools.@benchmark TS.mul!($res, $a, $b) setup=(TS.zero!($res)) evals=1 seconds=SECONDS samples=SAMPLES
    report(label, trial)
    return trial
end

function main()
    println("TaylorSeries product benchmarks")
    println("  DA order:    ", DA_ORDER)
    println("  time order:  ", TIME_ORDER)
    println("  seconds:     ", SECONDS)
    println("  samples:     ", SAMPLES)
    println()

    inputs = product_inputs()

    bench_product("TaylorN{Float64} product", inputs.f, inputs.g)
    bench_mul!("TaylorN{Float64} mul!", inputs.resN, inputs.f, inputs.g)
    bench_product("Taylor1{TaylorN{Float64}} product", inputs.tf, inputs.tg)
    bench_mul!("Taylor1{TaylorN{Float64}} mul!", inputs.resTn, inputs.tf, inputs.tg)
    bench_product("TaylorN{Taylor1{Float64}} product", inputs.fn_t1, inputs.gn_t1)
    bench_mul!("TaylorN{Taylor1{Float64}} mul!", inputs.resNt1, inputs.fn_t1, inputs.gn_t1)

    return nothing
end

main()
