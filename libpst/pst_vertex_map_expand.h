#ifndef pst_vertex_map_expand_H
#define pst_vertex_map_expand_H

/* pst_vertex_map_expand.h -- procedures for expanding vertex-located maps. */
/* Last edited on 2025-02-28 06:17:56 by stolfi */

#include <stdint.h>
#include <float_image.h>

float_image_t *pst_vertex_map_expand
  ( float_image_t *B,
    int32_t wch,
    int32_t NXA,
    int32_t NYA,
    double scale
  );
  /* The approximate inverse of {pst_vertex_map_shrink}. Given a height
    map {B}, produces a map {A} that is twice as big in the {X} (col) and {Y}
    (row) directions.  Then multiplies all data samples by {scale}.
    
    Assumes that the map elements are associated with the VERTICES of
    the integer grid. Namely, each value {A[c,X',Y']} is a property
    measured around the point of the {XY} plane with coordinates
    {(X',Y')}; and each value {B[c,X,Y]} will be associated with the
    point with coordinates {(X,Y)}. This is the case, for example, with
    height maps as defined in {pst_height_map.h}.
    
    In fractional index terms, point {(x,y)} of {A}'s domain corresponds
    to point {(x/2,y/2)} of {B}'s domain. Thus, for all indices {X,Y}
    {B[c,X,Y]} and {A[c,2*X,2*Y]} refer to the corresponding points in
    their respective domains.
    
    The output map {A} will have {NXA} columns and {NYA} rows. Thus, if
    {B} has {NXB} cols and {NYB} rows, one should typically specify
    {NXA=2*NXB-1+sx} and {NYA=2*NYB-1+sy} where {sx,sy} may be 0
    or 1.
    
    A sample {A[c,X',Y']} with even {X'} and {Y'} and any {c} is copied
    from {B[c,X'/2,Y'/2]}. For every channel {c} except {wch}, if {X'}
    and/or {Y'} is odd, the sample {A[c,X',Y']} is the weighted average
    of the two or four samples {B[c,X'/2+rx,Y'/2+ry]} where {rx} is 0 if
    {X'} is even and varies in {0,1} if {X'} is odd, and likewise for
    {ry} and {Y'}.
    
    The averaging will use the reliability weight of each
    sample{B[c,X,Y]}, assumed to be {B[wch,X,Y]} if channel {wch}
    exists; if it does not, the reliability weight is assumed to be 1.0.
    If a sample {B[c,X,Y]} does not exist in {B}, its reliability
    weight is assumed to be zero. Any sample {A[c,X',Y']} with zero
    weight is considered completely unknown and thus excluded from the
    averaging.
    
    Apart from the size reduction, every output sample value
    {A[c,X',Y']} will be multiplied bt {scale}, for every channel {c}
    except {wch}.
    
    If {wch} is a valid channel index in {B}, the output weight sample
    {A[wch,X',Y']} will be a combination of the weights {B[wch,X,Y]}
    of the samples that are used to compute {A[c,X',Y']}.
    For this combination, {B[wch,X,Y]} is interpreted as the
    reciprocal of the estimated variance of the noise in the sample
    value {B[c,X,Y]}.
    
    If the computed new weight {A[wch,X',Y']} is zero, or any of
    {A[c,X',Y']} comes out as {NAN} or overflows, then {A[c,X',Y']} will be
    set to zero if {c=wch} and to {NAN} otherwise. */

#endif
