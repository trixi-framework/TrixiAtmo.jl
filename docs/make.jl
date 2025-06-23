using Documenter
using DocumenterInterLinks
using TrixiAtmo

# Fix for https://github.com/trixi-framework/Trixi.jl/issues/668
if (get(ENV, "CI", nothing) != "true") &&
   (get(ENV, "TRIXI_DOC_DEFAULT_ENVIRONMENT", nothing) != "true")
    push!(LOAD_PATH, dirname(@__DIR__))
end

# Provide external links to the Trixi.jl docs (project root and inventory file)
links = InterLinks("Trixi" => ("https://trixi-framework.github.io/Trixi.jl/stable/",
                               "https://trixi-framework.github.io/Trixi.jl/stable/objects.inv"))

# Define module-wide setups such that the respective modules are available in doctests
DocMeta.setdocmeta!(TrixiAtmo, :DocTestSetup, :(using TrixiAtmo);
                    recursive = true)

# Copy some files from the repository root directory to the docs and modify them
# as necessary
# Based on: https://github.com/ranocha/SummationByPartsOperators.jl/blob/0206a74140d5c6eb9921ca5021cb7bf2da1a306d/docs/make.jl#L27-L41
open(joinpath(@__DIR__, "src", "code_of_conduct.md"), "w") do io
    # Point to source license file
    println(io,
            """
            ```@meta
            EditURL = "https://github.com/trixi-framework/Trixi.jl/blob/main/CODE_OF_CONDUCT.md"
            ```
            """)
    # Write the modified contents
    println(io, "# [Code of Conduct](@id code-of-conduct)")
    println(io, "")
    for line in eachline(joinpath(dirname(@__DIR__), "CODE_OF_CONDUCT.md"))
        line = replace(line,
                       "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref trixi_atmo_authors)")
        println(io, "> ", line)
    end
end

open(joinpath(@__DIR__, "src", "contributing.md"), "w") do io
    # Point to source license file
    println(io,
            """
            ```@meta
            EditURL = "https://github.com/trixi-framework/Trixi.jl/blob/main/CONTRIBUTING.md"
            ```
            """)
    # Write the modified contents
    for line in eachline(joinpath(dirname(@__DIR__), "CONTRIBUTING.md"))
        line = replace(line, "[LICENSE.md](LICENSE.md)" => "[License](@ref)")
        line = replace(line,
                       "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref trixi_atmo_authors)")
        println(io, line)
    end
end

# Since the authors header exists twice we create a unique identifier for the docs section.
authors_text = read(joinpath(dirname(@__DIR__), "AUTHORS.md"), String)
authors_text = replace(authors_text,
                       "[LICENSE.md](LICENSE.md)" => "[License](@ref)",
                       "# Authors" => "# [Authors](@id trixi_atmo_authors)")
write(joinpath(@__DIR__, "src", "authors.md"), authors_text)

open(joinpath(@__DIR__, "src", "license.md"), "w") do io
    write(io, "# License\n\n")
    for line in eachline(joinpath(dirname(@__DIR__), "LICENSE.md"))
        line = replace(line,
                       "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref trixi_atmo_authors)")
        write(io, "> " * line * "\n")
    end
end

# Copy contents from README to the starting page to not need to synchronize it manually
readme_text = read(joinpath(dirname(@__DIR__), "README.md"), String)
readme_text = replace(readme_text,
                      "[LICENSE.md](LICENSE.md)" => "[License](@ref)",
                      "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref trixi_atmo_authors)",
                      "<p" => "```@raw html\n<p",
                      "p>" => "p>\n```",
                      r"\[comment\].*\n" => "")    # remove comments
write(joinpath(@__DIR__, "src", "index.md"), readme_text)

makedocs(;
         modules = [TrixiAtmo],
         repo = Remotes.GitHub("trixi-framework", "TrixiAtmo.jl"),
         sitename = "TrixiAtmo.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://trixi-framework.github.io/TrixiAtmo.jl",
                                  edit_link = "main",
                                  assets = String[],),
         pages = ["Home" => "index.md",
             "Reference" => "reference.md",
             "Authors" => "authors.md",
             "Contributing" => "contributing.md",
             "Code of Conduct" => "code_of_conduct.md",
             "License" => "license.md"],
         plugins = [links],)

deploydocs(;
           repo = "github.com/trixi-framework/TrixiAtmo.jl",
           devbranch = "main",
           push_preview = true)
