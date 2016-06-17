using Documenter, TaylorSeries

makedocs(
    modules = Module[TaylorSeries],
    doctest = false,
)

deploydocs(
    deps = Deps.pip("pygments", "mkdocs", "mkdocs-cinder", "python-markdown-math"),
    repo   = "github.com/JuliaDiff/TaylorSeries.jl.git",
    osname = "linux"
)
