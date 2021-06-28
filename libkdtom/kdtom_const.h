/* An internal k-d-tree node that has a constant value over the whole domain. */
/* Last edited on 2021-06-28 03:33:59 by jstolfi */

#ifndef kdtom_const_H
#define kdtom_const_H

#define _GNU_SOURCE
#include <stdio.h>

#include <kdtom.h>
#include <ppv_array.h>
#include <bool.h>

typedef struct kdtom_const_t 
  { kdtom_t head;       /* General node parameters. */
    ppv_sample_t val;   /* Sample value. */
  } kdtom_const_t;
  /* A record {T} of this type describes a block {T.v} of voxels whose
    values are all equal to {T.val}.
    
    The field {head} must be the first field in the record. */

ppv_sample_t kdtom_const_get_sample(kdtom_const_t *T, ppv_index_t ix[]);
  /* Obtains the sample {T.v[ix]}, that is, {t.val}, if {ix} is a valid index vector for
    {T}; otherwise bombs out. */

kdtom_const_t *kdtom_const_make(ppv_dim_t d, ppv_nbits_t bps, ppv_size_t [], ppv_sample_t val);
  /* Creates a constant node {T} with the given parameters: {d} axes, {bps} bits per sample,
    and sample value {val}, with {size[k]} voxels along each axis {k} in {0..d-1}. */

size_t kdtom_const_node_bytesize(ppv_dim_t d);
  /* Size in bytes of a {kdtom_const_t} node {T}, including the {T.head.size}
    vector. */

size_t kdtom_const_bytesize(kdtom_const_t *T);
  /* Returns the size in bytes used by the record {*T}. The result
    includes the storage used by the {T.head} header and its
    {T.head.size} vector, as well as any internal type-specific
    fields. */

#endif
