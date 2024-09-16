"""
    MetricTermsCrossProduct()

Struct used for multiple dispatch on functions that compute the metric terms.
When the argument `metric_terms` is of type `MetricTermsCrossProduct`, the 
contravariant vectors are computed using the cross-product form.
"""
struct MetricTermsCrossProduct end

"""
    MetricTermsCurlInvariant()

Struct used for multiple dispatch on functions that compute the metric terms.
When the argument `metric_terms` is of type `MetricTermsCurlInvariant`, the 
contravariant vectors are computed using the curl invariant form.
"""
struct MetricTermsCurlInvariant end
