using Documenter, QSpin

makedocs(
    checkdocs = :exports,
    linkcheck = true,
    modules = [QSpin],
    pages = [
        "QSpin Home" => "index.md",
        "Overview" => "overview.md",
        "API Reference" => "api.md",
    ],
    remotes = nothing,
    sitename = "QSpin Documentation",
)
