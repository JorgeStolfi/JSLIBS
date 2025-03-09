#ifndef pst_vertex_map_shrink_H
#define pst_vertex_map_shrink_H

/* pst_vertex_map_shrink.h -- procedures for shrinking vertex-located maps. */
/* Last edited on 2025-02-27 10:21:27 by stolfi */

#include <stdint.h>
#include <float_image.h>

float_image_t *pst_vertex_map_shrink
  ( float_image_t *A,
    int32_t wch,
    int32_t NXB,
    int32_t NYB,
    double scale
  );
  /* Given a map {A}, returns another map {B} with same channel count,
    reduced by 1/2 in the {X} (col) and {Y} (row) directions. Then multiplies 
    all data samples in {B} by {scale}.
    
    Assumes that the map elements are associated with the VERTICES of
    the integer grid. Namely, each value {A[c,X',Y']} is a property
    measured around the point of the {XY} plane with coordinates
    {(X',Y')}; and, likewise, each value {B[c,X,Y]} will be associated with the
    point with coordinates {(X,Y)}. This is the case, for example, of
    height maps as defined in {pst_height_map.h}.
    
    The output map {B} will have {NXB} cols and {NYB} rows. If {A} has
    {NXA} cols and {NYA} rows, one should typically specify
    {NXB=NXA/2+1} and {NYB=NYA/2+1}, rounded DOWN.
    
    In fractional index terms, point {(x,y)} of {A}'s domain corresponds
    to point {(x/2,y/2)} of {B}'s domain. Thus, for all indices {X,Y}
    {B[c,X,Y]} and {A[c,2*X,2*Y]} refer to the corresponding points in
    their respective domains. Thus, for every channel {c} except {wch},
    sample {B[c,X,Y]} will be the weighted average of {A[c,X',Y']} where
    {X'=2*X+rx}, {Y'=2*Y+ry}, and {rx,ry} range in {-1..1}.
    
    The averaging will use the reliability weight of each sample
    {A[c,X',Y']}, assumed to {A[wch,X',Y']} if channel {wch} exists; if
    it does not, the reliability weight is assumed to be 1.0. If a
    sample {A[c,X',Y']} does not exist in {A}, its reliability weight is
    assumed to be zero. In either case this reliability weight is
    multiplied by {wH[rx]*wH[ry]} where {wH[0]=1} and {wH[±1]=1/2}. Any
    sample {A[c,X',Y']} with zero weight is considered completely
    unknown and thus excluded from the averaging.
    
    Apart from the size reduction, every output sample value {B[c,X,Y]}
    will be multiplied bt {scale}, for every channel {c} except {wch}.
    
    If {wch} is a valid channel index in {A}, the output weight sample
    {B[wch,X,Y]} will be a combination of the weights {A[wch,X',Y']}
    of the samples that are used to compute {B[c,X,Y]}.
    For this combination, {A[wch,X',Y']} is interpreted as the
    reciprocal of the estimated variance of the noise in the sample
    value {A[c,X',Y']}.
    
    If the computed new weight {B[wch,X,Y]} is zero, or any of
    {B[c,X,Y]} comes out as {NAN} or overflows, then {B[c,X,Y]} will be
    set to zero if {c=wch} and to {NAN} otherwise. */

#endif
