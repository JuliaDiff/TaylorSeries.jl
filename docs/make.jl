using Documenter, TaylorSeries

makedocs(
    modules = Module[TaylorSeries],
    clean   = false,
    doctest = false,
)

deploydocs(
    deps = Deps.pip("pygments", "mkdocs", "mkdocs-cinder", "python-markdown-math"),
    repo   = "github.com/JuliaDiff/TaylorSeries.jl.git",
)

