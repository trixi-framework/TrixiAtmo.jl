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

struct MetricTermsExactSpherical end

struct MetricTermsExactCartesian end