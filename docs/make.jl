using Documenter, SWMF

#push!(LOAD_PATH,"../src/")

makedocs(
    sitename="SWMF",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)

deploydocs(
    repo = "github.com/henry2004y/SWMF.git",
    target = "build",
    branch = "gh-pages"
)
