using Documenter
using DocumenterInterLinks
using TrixiAtmo

# Fix for https://github.com/trixi-framework/Trixi.jl/issues/668
if (get(ENV, "CI", nothing) != "true") &&
   (get(ENV, "TRIXI_DOC_DEFAULT_ENVIRONMENT", nothing) != "true")
    push!(LOAD_PATH, dirname(@__DIR__))
end

# Define module-wide setups such that the respective modules are available in doctests
DocMeta.setdocmeta!(TrixiAtmo, :DocTestSetup, :(using TrixiAtmo);
                    recursive = true)

# Provide external links to the Trixi.jl docs (project root and inventory file)
links = InterLinks("Trixi" => ("https://trixi-framework.github.io/Trixi.jl/stable/",
                               "https://trixi-framework.github.io/Trixi.jl/stable/objects.inv"))

# Copy list of authors to not need to synchronize it manually.
# Since the authors header exists twice we create a unique identifier for the docs section.
authors_text = read(joinpath(dirname(@__DIR__), "AUTHORS.md"), String)
authors_text = replace(authors_text,
                       "[LICENSE.md](LICENSE.md)" => "[License](@ref)",
                       "# Authors" => "# [Authors](@id trixi_atmo_authors)")
write(joinpath(@__DIR__, "src", "authors.md"), authors_text)

# Copy contributing information to not need to synchronize it manually
contributing_text = read(joinpath(dirname(@__DIR__), "CONTRIBUTING.md"), String)
contributing_text = replace(contributing_text,
                            "[LICENSE.md](LICENSE.md)" => "[License](@ref)",
                            "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref trixi_atmo_authors)")
write(joinpath(@__DIR__, "src", "contributing.md"), contributing_text)

# Copy code of conduct to not need to synchronize it manually
code_of_conduct_text = read(joinpath(dirname(@__DIR__), "CODE_OF_CONDUCT.md"), String)
code_of_conduct_text = replace(code_of_conduct_text,
                               "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref trixi_atmo_authors)")
write(joinpath(@__DIR__, "src", "code_of_conduct.md"), code_of_conduct_text)

# Copy license not need to synchronize it manually
open(joinpath(@__DIR__, "src", "license.md"), "w") do license_file
    write(license_file, "# License\n\n")
    for line in eachline(joinpath(dirname(@__DIR__), "LICENSE.md"))
        line_replaced = replace(line,
            "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref trixi_atmo_authors)")
        write(license_file, "> " * line_replaced * "\n")
    end
end

# Copy contents form README to the starting page to not need to synchronize it manually
readme_text = read(joinpath(dirname(@__DIR__), "README.md"), String)
readme_text = replace(readme_text,
                      "[LICENSE.md](LICENSE.md)" => "[License](@ref)",
                      "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref trixi_atmo_authors)",
                      "<p" => "```@raw html\n<p",
                      "p>" => "p>\n```",
                      r"\[comment\].*\n" => "")    # remove comments
write(joinpath(@__DIR__, "src", "home.md"), readme_text)

makedocs(;
         modules = [TrixiAtmo],
         authors = "Benedict Geihe <bgeihe@uni-koeln.de>, Tristan Montoya <montoya.tristan@gmail.com, Hendrik Ranocha <hendrik.ranocha@uni-mainz.de>, Andrés Rueda-Ramírez <am.rueda@upm.es>, Michael Schlottke-Lakemper <michael@sloede.com>",
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
