/* An inrernal k-d-tree node that divides the domain betwee two nodes. */
/* Last edited on 2024-12-05 10:32:58 by stolfi */

#ifndef kdtom_split_H
#define kdtom_split_H

#include <stdio.h>
#include <stdint.h>

#include <kdtom.h>
#include <ppv_array.h>
#include <bool.h>

typedef struct kdtom_split_t 
  { kdtom_t h;          /* General node parameters. */
    ppv_axis_t ax;      /* Coordinate axis perpendicular to cut. */
    ppv_size_t size0;     /* Size of low core domain part. */
    kdtom_t *sub[2];    /* Representation of the two half-blocks. */
  } kdtom_split_t;
  /* A record {T} of this type describes an infinite grid {T.V} of
    voxels whose core consists of the the splice of two grids {T0.V} and
    {T1.V}, separated by a plane perpendicular to axis {ax} between the
    voxels with indices {T.ixcut-1} and {T.ixcut} on that axis, where 
    {T.ixcut = T.h.ixlo[ax] + T->size0}. 
    
    The core {T.K} of {T} is defined by {T.h.ixlo} and {T.h.size}, as in
    any {kdtom_t} variant, and must not be empty. Any voxels outside
    that region are all equal to {T.fill}. Note that {T.fill},
    {T.sub[0].fill}, and {T.sub[1].fill} may be all different.

    The axis counts {T.sub[0].d} and {T.sub[1].d} must be equal, and
    that will be the count {T.h.d}. Same for the {.maxsmp} fields.
    
    The axis {T.ax} must be in {0..d-1} where {d = T.h.d}, and the value
    of {T.size0} must be in {0..T.d.size[ax]}. Note that the core {T.K}
    must not be empty (in particular, {T.h.size[ax]} must be at least 1).
    
    The core domain {T.DK} is split into two sub-domains:
    {T.DK0}, with with indices up to {T.ixcut-1} along axis {ax}, and
    {T.DK1} with indices from {T.ixcut} and higher. The value {T.K[ix]}
    of a voxel in {T}'s core is {T0.V[ix]} or {T1.V[ix]}, respectively.
    
    The sub-domain {T.DK0} will be empty if {size0 == 0}. Similarly
    {T.DK1} will be empty if {size0 == T.size[ax]} In either case the
    corresponding child pointer ({T.sub[0]} or {T.sub1}) will be {NULL}.
    This feature is useful if the other child has a different {fill}.
    Note that they cannot be both empty since {T.DK} is not empty.

    Any element {T0.V[ix]} with {ix} outside {T.DK0} is inaccessible
    through {T}. Thus normally the core domain {T0.DK} of {T0} is
    contained in the sub-domain {T.DK0}. Conversely, if {T0.DK} contains
    the whole of {T.DK0}, the fill value of {T0} is irrelevant. The same
    is true of {T1.V[ix]} relative to {T.DK1}.
     
    The grid {T0.V} is the grid {T.sub[0].V} implicitly translated by
    the low corner {T.ixlo0} of {T.DK0}. Similarly {T1.V} is {T.sub[1].V}
    implicitly translated by the low corner {T.ixlo1} of {T.DK1}. That is,
    a voxel {T.K[ix]} inside the core domain {T.DK} is either
    {T.sub[0].V[ix - T.ixlo0]} or {T.sub[1].V[ix - T.ixlo1]}, depending on
    whether {ix} is in {T.DK0} or {T.DK1}. The vector {T.ixlo0} is just
    {T.ixlo}, and {T.ixlo1} is the same except that {T.ixlo1[ax]} is
    {T.ixcut}.
    
    The {.d} and {.maxsmp} fields of {T.sub[0]} and {T.sub[1]} must be equal to 
    those of {T.d}.  The {.fill} values may be all different.
   
    The header field {h} must be the first field in the record. */

kdtom_split_t *kdtom_split_make
  ( ppv_dim_t d,
    ppv_sample_t maxsmp,
    ppv_sample_t fill,
    ppv_index_t ixlo[],
    ppv_size_t size[], 
    ppv_axis_t ax, 
    kdtom_t *sub0,
    ppv_size_t size0,
    kdtom_t *sub1,
    ppv_size_t size1
  );
  /* Creates a {kdtom_split_t} node {T} whose grid {T.V} is 
    the splice of the two grids {T0.V} and {T1.V}, clipped to the 
    sub-domains {T.DK0} and {T.DK1} defined by {ixlo}, {size},
    {size0}, and {size1}.
    
    This node will have {T.h.d, T.h.maxsmp, T.h.ixlo, T.h.size, T.size0}
    set as given.
    
    Exceptionally, if the domain {T.DK} is empty (that is, {size[k]==0}
    for any axis {k}), the result will be a {kdtom_const_t} node with
    empty core, everywhere equal to {fill}.
    
    If the domain {T.DK} is not empty, the node {T} will be a
    {kdtom_split_t} node with {T.sub[0] = sub0} and {T.sub[1] = sub1}.
    However, {sub0} must be {NULL} if and only if {size0==0}; and {sub1}
    must be {NULL} if and only if {size1==0}. They cannot be both
    {NULL}, because {size0+size1} must be {size[ax]}. The values of {d}
    and {maxsmp} must match those of the children nodes {sub0} and
    {sub1}, if they are not null.
    
    The procedure assumes that {sub0} and {sub1}, if not {NULL}, have
    been already clipped to the corresponding sub-domains ({T.DK0} and
    {T.DK1}) and have been translated so so that their core coordinates
    are relative to {T.ixlo0} and {T.ixlo1}, respectively.  
    
    Specifically, Assumes that, along axis {ax}, the indices in the core
    of {sub0.DK}, if any, are contained in {0..size0-1}, and those
    of{sub1.DK} are contained in 0..size1-1}. Along any other axis {k},
    the indices of the cores of either node, if any, must be in
    {0..size[k]-1}. */

kdtom_t *kdtom_split_make_clipping
  ( ppv_dim_t d,
    ppv_sample_t maxsmp,
    ppv_sample_t fill,
    ppv_index_t ixlo[], 
    ppv_size_t size[], 
    ppv_axis_t ax, 
    ppv_index_t ixcut,
    kdtom_t *T0,
    kdtom_t *T1
  );
  /* Tries to create a {kdtom_split_t} node {T} whose grid {T.V} is 
    the splice of the two grids {T0.V} and {T1.V}, clipped to the 
    sub-domains {T.DK0} and {T.DK1} defined by {ixlo}, {size},
    and {ixcut}.
    
    The values of {T.h.ixlo}, {T.h.size}, {T.h.fill} will be set as
    specified.  The two trees {T0} and {T1} must have the given
    values of {d} and {maxsmp}, which will be those of {T}.
    
    The typical case is when no element of {size} is zero, {size[ax]} is
    at least 2, and {ixcut-ixlo[ax]} is in {1..size[ax]-1}. Then in
    principle {T} will be a {kdtom_split_t} node, with {T.ax}, and
    {T.size0} set as described above. The fields {T.sub[0]} and {T.sub[1]}
    will be {T0} clipped to {T.DK0} and translated by {-T.ixlo0}, and {T1}
    clipped to {T.DK1} and translated by {-T.ixlo1}.
    
    Special cases
    
    However, for certain values of the parameters, this procedure may
    return a simpler {kdtom_t} node that is equivalent to the one
    described above.
    
    If any {size[k]} is zero, the core {T.K} will be empty, and {T.size}
    will be set to all zeros. Then the result {T} will be a
    {kdtom_const_t} node with the given {fill} value.  The parameters
    {T0,T1,ixcut} will be ignored.
    
    If {T.DK0} is empty (meaning that {ixcut <= ixlo[ax]}, so that
    {size0 == 0} and {size1 == size[ax]}), then {T0} is ignored (may be
    {NULL}), and the core {T.K} will be just {T1.V} clipped to {T.DK}
    and surrounded by {T1.fill}.
    
    Similarly, if {T.DK1} is empty (meaning {ixcut >=
    ixlo[ax]+size[ax]}, so that {size0 == size[ax]} and {size1 == 0}),
    then {T1} is ignored (may be {NULL}), and the core {T.K} will be
    just {T0.V} clipped to {T.DK} and surrounded by {T0.fill}. In both
    cases the value outside {T.DK} will be the given {fill} value.
    
    More generally, if the procedure determines that the grid {T.V} can
    be described by a single {kdtom_t} node with savings of space, the
    result will be such a node. */

kdtom_split_t *kdtom_split_clone(kdtom_split_t *T);
  /* Returns a copy of the node {T}. Copies the fixed fields, including the
    pointers {T->sub[0..1]}, and the internally allocated vectors {T.h.ixlo,T.h.size},
    but does NOT copy the nodes pointed to by {T->sub[0..1]}. */

kdtom_t *kdtom_split_simplify(kdtom_split_t *T);
  /* Tries to simplify the node {T}, for example by removing 
    {T.sub[0]} or {T.sub[1]} if they happen to be uniformly equal
    to {T.fill}, or if they can be joined into a single node.
    It may even replace the node {T} by one of those two nodes.
    
    The procedure will not modify {T}. If it does find a
    way to simplify it, returns a newly allocated node 
    with the changes. */

ppv_sample_t kdtom_split_get_core_sample(kdtom_split_t *T, ppv_index_t dx[]);
  /* Returns the sample {T.V[T.ixlo + dx]}. IMPORTANT: Assumes that this
    index is insde the core domain {T.DK}; that is, {dx} is in
    {0..T.size[k]-1}. The outcome is undefined otherwise. */

kdtom_t *kdtom_split_clip_core(kdtom_split_t *T, ppv_index_t ixlo[], ppv_size_t size[]);
  /* Returns a {kdtom_t} structure {S} that describes the same voxel
    grid as {T}, except that its core {S.K} will be the box {B} with
    low corner {ixlo[0..T.d-1]} and size {size[0..T.d-1]}.
    
    IMPORTANT: The procedures assumes that {B} is not empty and is
    contained in {T.DK}.
    
    The parameters {S.maxsmp} and {S.fill} will be copied from {T}.
    
    The result {S} will be a newly allocated node, even it describes the
    same grid as {T} or of its descendants. It may not have the same
    kind as {T}. The node {T} itself is not changed. Note that {T} and
    some of its descendants may not be reachable from {S}.
    
    The descendants of {T}, including the nodes {T.sub[0]} and {T.sub[1]}
    may have to be clipped too, recursively.  That may change the types
    of those descendants and/or cause whole sub-trees to be elminated. 
    
    In particular, if the box {B} does not intersect one of the core
    sub-domains {T.DK0} or {T.DK1}, the result will be just a copy of the other
    child node, respectively {T.sub[1]} or {T.sub[0]}, properly clipped 
    and translated. */

size_t kdtom_split_node_bytesize(ppv_dim_t d);
  /* Size in bytes of a {kdtom_split_t} node {T}, including the {T.h.ixlo} and {T.h.size}
    vector and the {T.sub[0..1]} pointers, but NOT the nodes in those
    sub-trees. The number of axes {d} should not be zero. */

size_t kdtom_split_bytesize(kdtom_split_t *T, bool_t total);
  /* If {total} is false, returns the size in bytes used by the node
    record {*T}. The result includes the storage used by the {T.h}
    header and the {T.h.ixlo} and {T.h.size} vectors, and the pointers
    {T.sub[0..1]}; but NOT the other nodes accessed through those two
    pointers.
    
    If {total} is true, the result also includes the total size in
    bytes of all nodes that can be reached through {T.sub[0..1]},
    as returned by {kdtom_bytesize(T.sub[i],TRUE)}. .*/

/* DEBUGGING */

void kdtom_split_print_fields(FILE *wr, kdtom_split_t *T);
  /* Prints to {wr} the fields that are specific of {T}, as ".{field} =
    {value}". */

void kdtom_split_print_subtrees(FILE *wr, int32_t ind, kdtom_split_t *T);
  /* Prints to {wr} the two sub-trees, with the roots indented
    by {ind} spaces. */

#endif
