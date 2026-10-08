# Extended benchmarks for the scalar-space / adopt changes (TaylorSeries.jl).
#
# Same usage as `benchmark_scalar_space.jl` (run on both versions, same machine):
#   julia --project -t1 benchmark_scalar_space_extended.jl save baseline.json
#   julia --project -t1 benchmark_scalar_space_extended.jl save patched.json
#   julia --project benchmark_scalar_space_extended.jl compare baseline.json patched.json [min]
# (the optional `min` compares minimum times instead of medians; the table also shows memory)
# Set BENCH_QUICK=1 for a faster, noisier run (0.5 s per benchmark instead of 2 s).

using BenchmarkTools, TaylorSeries

const QUICK = get(ENV, "BENCH_QUICK", "0") == "1"
const SECONDS = QUICK ? 0.5 : 2.0

# (TaylorN order, number of variables, Taylor1 order for the "light" Taylor1{TaylorN})
const CONFIGS = [
    (tag = "order6_2vars",  order = 6,  nv = 2, N = 10),
    (tag = "order10_3vars", order = 10, nv = 3, N = 10),
    (tag = "order6_6vars",  order = 6,  nv = 6, N = 10),
]

copy_t(a) = Taylor1([deepcopy(c) for c in a.coeffs])      # independent copy, for setup=

# Jet-transport workload: Picard iteration for Lorenz with Taylor1{TaylorN} variables
function picard_lorenz(x0, y0, z0, N)
    σ, ρ, β = 10.0, 28.0, 8/3
    x = Taylor1(x0, N); y = Taylor1(y0, N); z = Taylor1(z0, N)
    for _ in 1:N
        dx = σ * (y - x)
        dy = x * (ρ - z) - y
        dz = x * y - β * z
        x = x0 + integrate(dx)
        y = y0 + integrate(dy)
        z = z0 + integrate(dz)
    end
    return x, y, z
end

function make_suite(order, nv, N)
    variables!("x", order=order, numvars=nv, nowarn=true)
    xs = [TaylorN(Float64, i, order=order) for i in 1:nv]
    x, y = xs[1], xs[2]
    z = xs[min(3, nv)]
    s = sum(xs) / nv

    # light TaylorN (as in the basic file) and dense TaylorN (all monomials up to `order`)
    f  = 1 + x + 2y + x*y + 0.5x^2
    g  = 1 - y + 3x*y + y^2
    fd = exp(s)
    gd = cos(s) + x*y
    hp1 = f.coeffs[2]; hp2 = g.coeffs[2]
    hpd1 = fd.coeffs[order]; hpd2 = gd.coeffs[order]      # dense homogeneous polynomials

    tN  = Taylor1([f * i for i in 1:N+1], N)               # Taylor1{TaylorN}, light
    uN  = Taylor1([g * i for i in 1:N+1], N)
    Nd  = min(N, 5)                                         # dense Taylor1{TaylorN}: keep cost bounded
    tD  = Taylor1([fd * i for i in 1:Nd+1], Nd)
    uD  = Taylor1([gd * i for i in 1:Nd+1], Nd)
    t1  = Taylor1(N)                                        # Taylor1{Float64}
    cN  = convert(TaylorN{Float64}, 1.0)
    vN  = [f * i for i in 1:N+1]
    vals = fill(0.1, nv)

    suite = BenchmarkGroup()

    # --- 1. same-space operations, light and dense ----------------------------------
    s1 = suite["same-space ops"] = BenchmarkGroup()
    s1["light: TaylorN + TaylorN"]               = @benchmarkable $f + $g
    s1["light: TaylorN * TaylorN"]               = @benchmarkable $f * $g
    s1["light: HP + HP"]                         = @benchmarkable $hp1 + $hp1
    s1["light: HP * HP"]                         = @benchmarkable $hp1 * $hp2
    s1["light: Taylor1{TaylorN} + Taylor1{TaylorN}"] = @benchmarkable $tN + $uN
    s1["light: Taylor1{TaylorN} * Taylor1{TaylorN}"] = @benchmarkable $tN * $uN
    s1["light: zero(Taylor1{TaylorN})"]          = @benchmarkable zero($tN)
    s1["light: one(Taylor1{TaylorN})"]           = @benchmarkable one($tN)
    s1["dense: TaylorN + TaylorN"]               = @benchmarkable $fd + $gd
    s1["dense: TaylorN * TaylorN"]               = @benchmarkable $fd * $gd
    s1["dense: HP + HP"]                         = @benchmarkable $hpd1 + $hpd2
    s1["dense: HP * HP"]                         = @benchmarkable $hpd1 * $hpd2
    s1["dense: zero(TaylorN)"]                   = @benchmarkable zero($fd)
    s1["dense: Taylor1{TaylorN} + Taylor1{TaylorN}"] = @benchmarkable $tD + $uD
    s1["dense: Taylor1{TaylorN} * Taylor1{TaylorN}"] = @benchmarkable $tD * $uD

    # --- 2. numbers, converted constants, mixtures -------------------------------------
    s2 = suite["numbers and mixtures"] = BenchmarkGroup()
    s2["light: TaylorN + Float64"]               = @benchmarkable $f + 1.5
    s2["light: TaylorN + converted constant"]    = @benchmarkable $f + $cN
    s2["light: TaylorN * converted constant"]    = @benchmarkable $f * $cN
    s2["light: TaylorN{Int} + TaylorN{Float64}"] = @benchmarkable $(TaylorN(1, 2)) + $f
    s2["light: Taylor1{TaylorN} + Float64"]      = @benchmarkable $tN + 1.5
    s2["light: Taylor1{TaylorN} + Taylor1{Float64}"] = @benchmarkable $tN + $t1
    s2["light: Taylor1{TaylorN} * Taylor1{Float64}"] = @benchmarkable $tN * $t1
    s2["light: Taylor1{TaylorN} + TaylorN"]      = @benchmarkable $tN + $f
    s2["light: TaylorN * Taylor1{TaylorN}"]      = @benchmarkable $f * $tN
    s2["dense: TaylorN + Float64"]               = @benchmarkable $fd + 1.5
    s2["dense: TaylorN + converted constant"]    = @benchmarkable $fd + $cN
    s2["dense: Taylor1{TaylorN} + Taylor1{Float64}"] = @benchmarkable $tD + $(Taylor1(Nd))
    s2["dense: Taylor1{TaylorN} * TaylorN"]      = @benchmarkable $tD * $fd
    s2["convert(TaylorN{Float64}, 1.0)"]         = @benchmarkable convert(TaylorN{Float64}, 1.0)

    # --- 3. containers, copies, setindex! ----------------------------------------------
    s3 = suite["containers and setindex!"] = BenchmarkGroup()
    s3["[x, 1.0]"]                               = @benchmarkable [$x, 1.0]
    s3["[x, 1.0, ..., 1.0] (10 numbers)"]        = @benchmarkable [$x, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    s3["Taylor1(Vector{Float64}) (N=11)"]        = @benchmarkable Taylor1($(rand(N+1)))
    s3["Taylor1(Vector{TaylorN}) light (N+1)"]   = @benchmarkable Taylor1($vN)
    s3["Taylor1(Vector{TaylorN}, order) light"]  = @benchmarkable Taylor1($vN, $(N+3))
    s3["Taylor1([x, 1.0, 2.0])"]                 = @benchmarkable Taylor1([$x, 1.0, 2.0])
    s3["Taylor1(fill(f, N+1)) (repeated object)"] = @benchmarkable Taylor1(fill($f, $(N+1)))
    s3["Taylor1(fill(fd, N+1)) dense (repeated)"] = @benchmarkable Taylor1(fill($fd, $(N+1)))
    s3["Taylor1{TaylorN} t[k] = TaylorN"]        = @benchmarkable (a[3] = $g) setup=(a = copy_t($tN)) evals=1
    s3["Taylor1{TaylorN} N+1 stores in a loop"]  = @benchmarkable (for k in 0:$N; a[k] = $f; end) setup=(a = copy_t($tN)) evals=1
    s3["TaylorN a[k] = HP (light)"]              = @benchmarkable (a[1] = $hp1) setup=(a = deepcopy($f)) evals=1
    s3["TaylorN a[k] = HP (dense)"]              = @benchmarkable (a[$(order-1)] = $hpd1) setup=(a = deepcopy($fd)) evals=1
    s3["TaylorN(Vector{HP})"]                    = @benchmarkable TaylorN($(f.coeffs[:]))
    s3["TaylorN(Vector{HP}) dense"]              = @benchmarkable TaylorN($(fd.coeffs[:]))
    s3["convert(Taylor1{TaylorN{Float64}}, same type)"] = @benchmarkable convert(Taylor1{TaylorN{Float64}}, $tN)

    # --- 4. scaling with the number of coefficients (distinct objects) -----------------
    s4 = suite["scaling of Taylor1(Vector{TaylorN})"] = BenchmarkGroup()
    for n in (10, 30, 60, 100, 200)       # n > 64 uses an IdSet in the duplicate scan
        vn = [f * i for i in 1:n+1]
        s4["Taylor1(vector), distinct, n=$n"] = @benchmarkable Taylor1($vn)
        s4["Taylor1(fill), repeated,  n=$n"]  = @benchmarkable Taylor1(fill($f, $(n+1)))
    end

    # --- 5. functions, calculus, evaluation on Taylor1{TaylorN} -----------------------
    s5 = suite["functions, calculus, evaluation"] = BenchmarkGroup()
    s5["exp(Taylor1{TaylorN}) dense"]            = @benchmarkable exp($tD)
    s5["sin(Taylor1{TaylorN}) dense"]            = @benchmarkable sin($tD)
    s5["exp(TaylorN) dense"]                     = @benchmarkable exp($fd)
    s5["sqrt(1 + TaylorN) dense"]                = @benchmarkable sqrt(1 + 0.1 * $fd)
    s5["differentiate(Taylor1{TaylorN}) light"]  = @benchmarkable differentiate($tN)
    s5["integrate(Taylor1{TaylorN}) light"]      = @benchmarkable integrate($tN)
    s5["gradient(TaylorN) dense"]                = @benchmarkable TS.gradient($fd)
    s5["evaluate(TaylorN, point) dense"]         = @benchmarkable evaluate($fd, $vals)
    s5["evaluate(Taylor1{TaylorN}, number)"]     = @benchmarkable evaluate($tN, 0.1)
    s5["evaluate(Taylor1{TaylorN}, point)"]      = @benchmarkable evaluate($tN, $vals)

    # --- 6. workload: Picard iteration for Lorenz (jet transport) -----------------------
    s6 = suite["workload: Lorenz Picard iteration"] = BenchmarkGroup()
    x0 = 1.0 + xs[1]; y0 = 2.0 + xs[min(2, nv)]; z0 = 3.0 + z
    s6["Picard, N=4"]                            = @benchmarkable picard_lorenz($x0, $y0, $z0, 4)
    s6["Picard, N=6"]                            = @benchmarkable picard_lorenz($x0, $y0, $z0, 6)

    # --- 7. explicit (non-default) space: new behaviour --------------------------------
    s7 = suite["explicit space (new behaviour)"] = BenchmarkGroup()
    try
        sp = JetSpace(order=order, variables=["u$i" for i in 1:nv])
        us = variables(sp)
        u = us[1]
        fu = 1 + u + 2us[2] + u*us[2]
        tu = Taylor1([fu * i for i in 1:N+1], N)
        cu = convert(TaylorN{Float64}, 1.0)
        s7["TaylorN + converted constant"]        = @benchmarkable $fu + $cu
        s7["TaylorN * converted constant"]        = @benchmarkable $fu * $cu
        s7["[u, 1.0] -> Taylor1"]                 = @benchmarkable Taylor1([$u, 1.0])
        s7["Taylor1{TaylorN} + Taylor1{Float64}"] = @benchmarkable $tu + $t1
        s7["Taylor1(fill(u, N+1))"]               = @benchmarkable Taylor1(fill($u, $(N+1)))
    catch err
        @warn "explicit-space group skipped" err
    end
    return suite
end

# --- driver ---------------------------------------------------------------------------
function run_group!(results, gname, grp)
    haskey(results, gname) || (results[gname] = BenchmarkGroup())
    res = results[gname]
    for (bname, b) in grp
        try
            b.params.seconds = SECONDS
            tune!(b)
            res[bname] = run(b)
        catch err
            msg = first(sprint(showerror, err), 120)
            @warn "n/a: $gname / $bname" exception=msg
        end
    end
end

function run_all()
    results = BenchmarkGroup()
    for cfg in CONFIGS
        @info "configuration $(cfg.tag)" cfg
        cres = results[cfg.tag] = BenchmarkGroup()
        suite = try
            make_suite(cfg.order, cfg.nv, cfg.N)
        catch err
            @warn "configuration skipped: $(cfg.tag)" exception=first(sprint(showerror, err), 200)
            continue
        end
        for (gname, grp) in suite
            run_group!(cres, gname, grp)
        end
    end
    return results
end

function print_results(res)
    for tag in sort(collect(keys(res)))
        println("\n########## ", tag)
        for gname in sort(collect(keys(res[tag])))
            println("\n## ", gname)
            for bname in sort(collect(keys(res[tag][gname])))
                r = res[tag][gname][bname]
                println(rpad(bname, 62), rpad(BenchmarkTools.prettytime(time(median(r))), 12),
                    "allocs=", allocs(r), "   mem=", BenchmarkTools.prettymemory(memory(r)))
            end
        end
    end
end

# `est` is the estimator used for the times: `median` (default) or `minimum`, which is less
# sensitive to GC pauses and machine noise for deterministic code (try it on the large
# 6-variable cases, where each call allocates tens of kB).
function compare(old, new; est = median)
    println(rpad("benchmark", 64), rpad("old", 12), rpad("new", 12), rpad("ratio", 9),
        rpad("allocs old->new", 18), "memory old->new")
    for tag in sort(collect(keys(new)))
        println("\n########## ", tag)
        for gname in sort(collect(keys(new[tag])))
            println("\n## ", gname)
            for bname in sort(collect(keys(new[tag][gname])))
                n = new[tag][gname][bname]
                if haskey(old, tag) && haskey(old[tag], gname) && haskey(old[tag][gname], bname)
                    o = old[tag][gname][bname]
                    ratio = time(est(n)) / time(est(o))
                    println(rpad(bname, 64), rpad(BenchmarkTools.prettytime(time(est(o))), 12),
                        rpad(BenchmarkTools.prettytime(time(est(n))), 12),
                        rpad(string(round(ratio, digits=2), "x"), 9),
                        rpad(string(allocs(o), " -> ", allocs(n)), 18),
                        BenchmarkTools.prettymemory(memory(o)), " -> ", BenchmarkTools.prettymemory(memory(n)))
                else
                    println(rpad(bname, 64), rpad("n/a", 12),
                        rpad(BenchmarkTools.prettytime(time(est(n))), 12), "(new only)")
                end
            end
        end
    end
    println("\nNoise is typically 5-10%; look at ratios clearly above ~1.15 and at allocation counts.")
end

function main(args)
    isempty(args) && (args = ["run"])
    cmd = args[1]
    if cmd == "run" || cmd == "save"
        res = run_all()
        print_results(res)
        cmd == "save" && (BenchmarkTools.save(args[2], res); println("\nsaved ", args[2]))
    elseif cmd == "compare"
        est = length(args) >= 4 && args[4] == "min" ? minimum : median
        compare(BenchmarkTools.load(args[2])[1], BenchmarkTools.load(args[3])[1]; est)
    else
        println("usage: run | save FILE | compare OLD NEW")
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    try
        main(ARGS)
    catch err
        Base.display_error(err, catch_backtrace())   # the real error, not just this call
        exit(1)
    end
end