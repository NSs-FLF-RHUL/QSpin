using Documenter, QSpin

makedocs(
    sitename = "QSpin Documentation",
    remotes = nothing,
    modules = [QSpin],
    checkdocs = :exports,
)
