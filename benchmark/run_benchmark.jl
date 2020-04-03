using PkgBenchmark
results = benchmarkpkg("TaylorSeries")
show(results)

#=
# specify tag and uncommit to benchmark versus prior tagged version
tag =
results = judge("TaylorSeries", tag)
show(results)
=#

export_markdown("..\\results.md", results)
