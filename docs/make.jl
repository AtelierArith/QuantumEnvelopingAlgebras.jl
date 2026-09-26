using Documenter
using QuantumEnvelopingAlgebras

makedocs(;
    sitename = "QuantumEnvelopingAlgebras.jl",
    modules = [QuantumEnvelopingAlgebras],
    authors = "Satoshi Terasaki",
    remotes = nothing,
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "API reference" => "api.md",
    ],
)

# Uncomment and set the repository URL once it is hosted, then remove
# `remotes = nothing` above so source links are generated.
#
# deploydocs(; repo = "github.com/USER_NAME/QuantumEnvelopingAlgebras.jl.git")
