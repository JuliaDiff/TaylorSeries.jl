# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

using PrecompileTools

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    t = Taylor1(20)
    δ = set_variables("δ", order=6, numvars=2)
    tN = one(δ[1]) + Taylor1(typeof(δ[1]), 20)
    # tb = Taylor1(Float128, 20)
    # δb = zero(Float128) .+ δ
    # tbN = one(δb[1]) + Taylor1(typeof(δb[1]), 20)
    #
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        for x in (:t, :tN)
            @eval begin
                T = numtype($x)
                zero($x)
                sin($x)
                cos($x)
                $x/sqrt($x^2+(2*$x)^2)
                evaluate(($x)^3, 0.125)
                ($x)[2]
            end
        end
    end
end
