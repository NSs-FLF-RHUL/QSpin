using Documenter, QSpin

makedocs(
    sitename = "QSpin Documentation",
    remotes = nothing,
    modules = [QSpin],
    linkcheck = true,
    checkdocs = :exports,
)
