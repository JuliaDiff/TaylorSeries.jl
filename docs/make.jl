using Documenter, TaylorSeries

makedocs(
    modules = [TaylorSeries],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "TaylorSeries.jl",
    authors  = "Luis Benet and David P. Sanders",
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
    deps   = nothing,
    make   = nothing,
    push_preview = true
)
