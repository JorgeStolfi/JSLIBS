/* Converting a {ppv_array_t} into a parsimonious k-d-tree. */
/* Last edited on 2024-12-05 10:32:53 by stolfi */

#ifndef kdtom_grind_array_H
#define kdtom_grind_array_H

#include <stdio.h>

#include <ppv_array.h>
#include <bool.h>

#include <kdtom.h>

kdtom_t *kdtom_grind_array(ppv_array_t *A, ppv_sample_t fill);
  /* Creates a {kdtom_t} structure {T} whose core {T.K} 
    is a k-d-tree derived from the array {A} by recursively splitting
    and analyzing each part.
    
    The non-leaf nodes of the tree, if any, will be {kdtom_split_t} nodes. The
    leaves will be either {kdtom_const_t} or {kdtom_array_t} nodes.
    The parameter {T.maxsmp} will be {A.maxsmp}, and {T.fill} will
    be the given {fill} value.
    
    Any substantial part of {A} found in the search whose voxels have
    all the same value is represented by a {kdtom_const_t} node.
    Otherwise that part of {A} is represented by a {kdtom_array_t}
    that shares storage with {A}.
    
    The recursive subdivision will proceed until the part of {A} under
    analysis would take more space if represented as a k-d-tree than if
    it was kept as a single array node.
    
    In particular, if all elements of {A} have the same value, the
    result is a single {kdtom_const_t} node with that value. If the two
    sides of a split turn out to be {kdtom_const_t} nodes with compatible
    parameters they are merged into one. Ditto if they are both
    {kdtom_array_t} nodes.
    
    IMPORTANT: since the storage areas of the {kdtom_array_t} nodes are
    shared parts of {A}'s storage area, neither they nor {A} can be
    wholly reclaimed. To elminate such sharing, use
    {kdtom_realloc_array_nodes} below. */

#endif

