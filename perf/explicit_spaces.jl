# Focused checks for explicit JetSpaces.
#
# Run from the repository root with:
#     julia --project=. perf/explicit_spaces.jl

using TaylorSeries

function report(label, f, args...)
    f(args...)
    GC.gc()
    alloc = @allocated f(args...)
    GC.gc()
    elapsed = @elapsed f(args...)
    println(rpad(label, 34), " time=", elapsed, "s alloc=", alloc, " bytes")
    return nothing
end

bench_multiply(f, g) = f * g

function bench_mul_loop!(res, f, g, n)
    for _ in 1:n
        TS.zero!(res)
        for k in eachindex(res)
            TS.mul!(res, f, g, k)
        end
    end
    return res
end

space = JetSpace(order=6, variables=[:x1, :x2, :x3, :x4])
x = get_variables(space)

f = sin(x[1] + x[2]*x[3]) + exp(x[4])
g = cos(x[2] - x[4]) + x[1]^2
res = zero(f)

report("TaylorN multiplication", bench_multiply, f, g)

report("Repeated mul! loop", bench_mul_loop!, res, f, g, 100)

du = [zero(f) for _ in 1:4]
u = [1.0 + x[1], -0.5 + x[2], 0.25 + x[3], x[4]]

function rhs!(du, u)
    du[1] = u[2] + sin(u[1]*u[3])
    du[2] = -u[1] + u[4]^2
    du[3] = u[1]*u[2] - cos(u[4])
    du[4] = exp(u[1] - u[3])
    return du
end

report("RHS-style vector field", rhs!, du, u)

orderT = 8
tu = [Taylor1([ui], orderT) for ui in u]
tdu = [zero(tui) for tui in tu]

function jet_rhs!(du, u)
    du[1] = u[2] + sin(u[1]*u[3])
    du[2] = -u[1] + u[4]^2
    du[3] = u[1]*u[2] - cos(u[4])
    du[4] = exp(u[1] - u[3])
    return du
end

report("Small Taylor1{TaylorN} RHS", jet_rhs!, tdu, tu)
