using Documenter
using DocumenterInterLinks

# Fix for https://github.com/trixi-framework/Trixi.jl/issues/668
if (get(ENV, "CI", nothing) != "true") &&
   (get(ENV, "TRIXI_DOC_DEFAULT_ENVIRONMENT", nothing) != "true")
    push!(LOAD_PATH, dirname(@__DIR__))
end

using TrixiAtmo

# Provide external links to the Trixi.jl docs (project root and inventory file)
links = InterLinks("Trixi" => ("https://trixi-framework.github.io/Trixi.jl/stable/",
                               "https://trixi-framework.github.io/Trixi.jl/stable/objects.inv"))

DocMeta.setdocmeta!(TrixiAtmo, :DocTestSetup, :(using TrixiAtmo);
                    recursive = true)

makedocs(;
         modules = [TrixiAtmo],
         authors = "Benedict Geihe <bgeihe@uni-koeln.de>, Tristan Montoya <montoya.tristan@gmail.com, Hendrik Ranocha <hendrik.ranocha@uni-mainz.de>, Michael Schlottke-Lakemper <michael@sloede.com>",
         repo = Remotes.GitHub("trixi-framework",
                               "TrixiAtmo.jl/blob/{commit}{path}#{line}"),
         sitename = "TrixiAtmo.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://trixi-framework.github.io/TrixiAtmo.jl",
                                  edit_link = "main",
                                  assets = String[],),
         pages = ["Home" => "index.md"],
         plugins = [links],)

deploydocs(;
           repo = "github.com/trixi-framework/TrixiAtmo.jl",
           devbranch = "main",
           push_preview = true)
