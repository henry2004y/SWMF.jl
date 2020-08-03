using Documenter, Batsrus

#push!(LOAD_PATH,"../src/")

makedocs(
    sitename="Batsrus",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)

deploydocs(
    repo = "github.com/henry2004y/Batsrus.jl.git",
    target = "build",
    branch = "gh-pages"
)
