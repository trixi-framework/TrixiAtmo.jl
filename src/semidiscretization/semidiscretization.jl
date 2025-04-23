"""
    MetricTermsCrossProduct()

Struct used for multiple dispatch on functions that compute the metric terms.
When the argument `metric_terms` is of type `MetricTermsCrossProduct`, the 
contravariant vectors are computed using the cross-product form.
"""
struct MetricTermsCrossProduct end

"""
    MetricTermsInvariantCurl()

Struct used for multiple dispatch on functions that compute the metric terms.
When the argument `metric_terms` is of type `MetricTermsInvariantCurl`, the 
contravariant vectors are computed using the invariant curl form.
## References
- Kopriva, D. A. (2006). Metric identities and the discontinuous spectral element method on 
  curvilinear meshes. Journal of Scientific Computing 26, 301-327. 
  [DOI: 10.1007/s10915-005-9070-8](https://doi.org/10.1007/s10915-005-9070-8)
- Vinokur, M. and Yee, H. C. (2001). Extension of efficient low dissipation high order
  schemes for 3-D curvilinear moving grids. In Caughey, D. A., and Hafez, M. M. (eds.),
  Frontiers of Computational Fluid Dynamics 2002, World Scientific, Singapore, pp. 129–164.
  [DOI: 10.1142/9789812810793_0008](https://doi.org/10.1142/9789812810793_0008)
"""
struct MetricTermsInvariantCurl end

"""
    MetricTermsExactSphere(ChristoffelSymbolsAutodiff())

Struct specifying options for computing geometric information for discretizations in 
covariant form based on an exact representation of the spherical geometry.  Currently, the 
only field is `christoffel_symbols`, specifying the approach used to compute the Christoffel symbols, for which the options are [`ChristoffelSymbolsAutodiff`](@ref) or 
[`ChristoffelSymbolsCollocationDerivative`](@ref). 
"""
struct MetricTermsCovariantSphere{ChristoffelSymbols}
    christoffel_symbols::ChristoffelSymbols
    function MetricTermsCovariantSphere(christoffel_symbols::ChristoffelSymbols = ChristoffelSymbolsAutodiff()) where {ChristoffelSymbols}
        return new{ChristoffelSymbols}(christoffel_symbols)
    end
end

@doc raw"""
    ChristoffelSymbolsAutodiff()

Struct used for multiple dispatch on functions that compute the Christoffel symbols. This 
option uses [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) to compute
```math
\Gamma_{jk}^i =
\frac{1}{2}G^{il}\big(\partial_j G_{kl} + \partial_k G_{jl} - \partial_l G_{jk}\big)
```
using forward-mode automatic differentiation. 
!!! warning
    Using this option with [`GlobalSphericalCoordinates`](@ref) is prone to `NaN` values as 
    a result of the polar singularity. This is remedied through the use of 
    [`GlobalCartesianCoordinates`](@ref).
"""
struct ChristoffelSymbolsAutodiff end

@doc raw"""
    ChristoffelSymbolsCollocationDerivative()

Struct used for multiple dispatch on functions that compute the Christoffel symbols. 
Letting $I^N$ denotes a degree $N$ polynomial interpolation operator on the scheme's quadrature nodes, this option computes the Christoffel symbols at each quadrature node 
using the approximation
```math
\Gamma_{jk}^i \approx
\frac{1}{2}G^{il}\big(\partial_j I^N G_{kl} + \partial_k I^N G_{jl} - \partial_l I^N G_{jk}\big).
```
"""
struct ChristoffelSymbolsCollocationDerivative end

include("semidiscretization_hyperbolic_2d_manifold_in_3d.jl")
