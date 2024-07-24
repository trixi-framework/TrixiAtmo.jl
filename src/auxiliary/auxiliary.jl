"""
    examples_dir()

Return the directory where the example files provided with TrixiAtmo.jl are located.
If TrixiAtmo.jl is installed as a regular package (with `]add Trixi`), these files are
read-only and should *not* be modified.
To find out which files are available, use, e.g., `readdir`:

# Examples
```@example
readdir(examples_dir())
```
"""
examples_dir() = pkgdir(TrixiAtmo, "examples")
