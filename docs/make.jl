using Documenter
using DataEnabledControl

makedocs(
    sitename = "DataEnabledControl",
    authors= "Mohammad <msramada@eng.ucsd.edu>",
    format = Documenter.HTML(),
    modules = [DataEnabledControl],
    pages=[
        "Home" => "index.md",
        "Theory" => "theory.md",
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/msramada/DataEnabledControl"
)
