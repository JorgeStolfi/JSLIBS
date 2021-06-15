/* Brushes (structuring elements) represented as lists of voxel indices. */
/* Last edited on 2021-06-13 13:36:14 by jstolfi */

#ifndef ppv_brush_H
#define ppv_brush_H

#define _GNU_SOURCE
#include <stdint.h>

#include <ppv_array.h>

typedef struct ppv_brush_t ppv_brush_t;
  /* A data structure describing a /brush/ -- a set {V(b)} of voxel
    index tuples relative to some central voxel. A brush is meant to
    define the (unweighted) neighborhood set of a position-invariant
    local operator, like erosion dilation, median filtering, etc. It is
    a "structuring element" in the mathematical morphlogy nomenclature.
    
    The central voxel, by definition, has relative indices {[0,0,...]}
    and thus the indices in {V(b)} may be negative.  The central
    voxel may not be in {V(b)}. */

ppv_dim_t ppv_brush_dimension(ppv_brush_t *b);
  /* Returns the dimension (number of indices, number of axes) of brush {b}. */

ppv_sample_count_t ppv_brush_voxel_count(ppv_brush_t *b);
  /* Returns the total count of voxels in the brush {b}. */
  
void ppv_brush_index_ranges(ppv_brush_t *b, ppv_index_t ixlo[], ppv_index_t ixhi[]);
  /* For each axis {ax} in {0..d-1}, sets {ixlo[ax]} and {ixhi[ax]} to
    the min and max values of the voxel indices of {b} along that axis;
    where {d} is {b}'s dimension. */

bool_t ppv_brush_enum(ppv_index_op_t op, ppv_brush_t *b);
  /* Calls {op(ix)} for all index tuples {ix} in the brush {b}.
    If any calll to {op} returns {TRUE}, stops the enumeration 
    immediately and returns {TRUE}.  Otherwise returns {FALSE}. */

void ppv_brush_free(ppv_brush_t *b);
  /* Reclaims the storage used by the internal tables of {b}, including {b} itself. */

/* CREATING BRUSHES */

#define ppv_brush_ball_MAX_RADIUS (100.0)
  /* An arbitrary limit on the radius of a ball brush. */

ppv_brush_t *ppv_brush_make_ball(ppv_dim_t d, double rad);
  /* Returns a {ppv_brush_t} structure for a brush that is a
    digital ball of radius {rad}. The radius must be non-negative
    and must not exceed {ppv_brush_ball_MAX_RADIUS}.
    
    The brush set {V(b)} will be the indices of all voxels of the
    bi-infinte integer voxel grid whose centers lie at distance {|rad|}
    (in voxels) or less from the voxel with indices {{0,0,...}}. The
    resulting brush will be symmetric about the center voxel.
    
    If {rad} is zero, the ball will have a single voxel, the central
    one. In general, the number of voxels will be approximately
    {4.19*|rad|^3}. */

#endif
