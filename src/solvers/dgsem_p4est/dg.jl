# Extract contravariant vector Ja^i (i = index) as SVector
# This function dispatches on the type of contravariant_vectors
static2val(::Trixi.StaticInt{N}) where {N} = Val{N}()
@inline function Trixi.get_contravariant_vector(index, contravariant_vectors::PtrArray,
                                                indices...)
    SVector(ntuple(@inline(dim->contravariant_vectors[dim, index, indices...]),
                   static2val(static_size(contravariant_vectors, Trixi.StaticInt(1)))))
end
