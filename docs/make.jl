using Documenter
using DataEnabledControl

makedocs(
    sitename = "DataEnabledControl",
    format = Documenter.HTML(),
    modules = [DataEnabledControl]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
