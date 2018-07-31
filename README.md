# Tensor Transforms

An arbitrary matrix can be used to transform n dimensional objects, however only a subset of these transformations
can be applied to a tensor in a way that makes sense. Specifically, any transform that can be represented by a signed
permutation matrix can be applied to a tensor in this way, and the set of signed permutation matrices is sufficient to
represent any and all possible rotations, inversions, and transposes of tensors.

This library provides functions for applying, manipulating, and generating this type of transformation. In the future
it will provide optional methods to apply these transformations on objects from other rust math libraries.