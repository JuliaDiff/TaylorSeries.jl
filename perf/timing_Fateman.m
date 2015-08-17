(* Timing a test by Fateman

In a mac, run this as
/Applications/Mathematica.app/Contents/MacOS/MathematicaScript -script ./timing_Fateman.m)

The `timeit` function is slightly modified from
https://github.com/JuliaLang/julia/blob/master/test/perf/micro/perf.nb
which is licensed under MIT

*)

ClearAll[timeit];
SetAttributes[timeit, HoldFirst];
timeit[ex_, name_String] := Module[
        {t},
        t = Infinity;
        Do[
                t = Min[t, N[First[AbsoluteTiming[ex]]]];
                ,
                {i, 1, 5}
        ];
    Print["mathematica,", name, ",min(time),", t];
];

(* Print[First[AbsoluteTiming[Function[n,Expand[(1+x+y+z+w)^n * (1+(1+x+y+z+w)^n)]][20]]]] *)

timeit[Function[n,Expand[(1+x+y+z+w)^n * (1+(1+x+y+z+w)^n)]][20],"fateman"]
