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

@doc raw"""
    cartesian_to_sphere(x)

Convert the global Cartesian coordinates `x` of a point in three-dimensional space into
spherical coordinates, returning the tuple `(lambda, phi, r)` of longitude
``\lambda \in [0, 2\pi)``, latitude ``\varphi \in [-\pi/2, \pi/2]``, and distance ``r`` from
the origin.
"""
@inline function cartesian_to_sphere(x)
    r = norm(x)
    lambda = atan(x[2], x[1])
    if lambda < 0
        lambda += 2 * pi
    end
    phi = asin(x[3] / r)

    return lambda, phi, r
end
