using Documenter, TaylorSeries

makedocs(
    # modules = TaylorSeries,
    clean   = false,
    doctest = false,
)

# deploydocs(
#     deps = Deps.pip("pygments", "mkdocs", "mkdocs-bootstrap", "python-markdown-math"),
#     repo   = "github.com/JuliaDiff/TaylorSeries.jl.git",
#     julia  = "0.4",
#     osname = "osx"
# )
