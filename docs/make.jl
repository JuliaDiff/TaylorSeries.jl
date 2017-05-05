using Documenter, TaylorSeries

makedocs(
    modules = [TaylorSeries],
    format = :html,
    sitename = "TaylorSeries.jl",
    pages = [
        "Home" => "index.md",
        "Background" => "background.md",
        "User guide" => "userguide.md",
        "Examples" => "examples.md",
        "API" => "api.md"
    ]
)

deploydocs(
    repo   = "github.com/JuliaDiff/TaylorSeries.jl.git",
    target = "build",
    julia = "0.5",
    osname = "linux",
    deps   = nothing,
    make   = nothing
)
