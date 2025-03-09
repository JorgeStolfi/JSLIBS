#ifndef pst_cell_map_shrink_H
#define pst_cell_map_shrink_H

/* pst_cell_map_shrink.h -- procedures for working with slope maps. */
/* Last edited on 2025-02-27 14:13:13 by stolfi */

#include <stdint.h>
#include <bool.h>

#include <float_image.h>

float_image_t *pst_cell_map_shrink
  ( float_image_t *A,
    int32_t wch,
    int32_t NXB,
    int32_t NYB,
    double scale
  );
  /* Given a map {A}, returns another map {JG},
    with half the size as {A} in the {X} (col) and {Y} (row)
    directions.  Then multiplies all data samples in {B} by {scale}.
    
    Assumes that the map elements are associated with the CELLS of the
    integer grid. Namely, each value {A[c,X',Y']} is an average of some
    quantity over the unit square of the {XY} plane with lower corner
    {(X',Y')} and center {(X'+1/2,Y'+1/2)}; and, likewise, each value
    {B[c,X,Y]} will be the averge over the unit square with lower corner
    {(X,Y)} and center {(X+1/2,Y+1/2)}. This is the case, for example,
    of slope maps as defined in {pst_slope_map.h}.
     
    The output map {B} will have {NXB} cols and {NYB} rows. If {A} has
    {NXA} cols and {NYA} rows, one should typically specify {NXB=NXA/2}
    and {NYB=NYA/2}, rounded UP.
   
    In fractional index terms, point {(x,y)} of {A}'s domain corresponds
    to point {(x/2,y/2)} of {B}'s domain. Thus, for all indices {X,Y},
    the cell of sample {B[c,X,Y]} in {B}'s domain corresponds to the
    {2×2} square in {A}'s domain with lower corner {(2*X,2*Y)}. Thus,
    for every channel {c} except {wch}, sample {B[c,X,Y]} will be the
    weighted average of {A[c,X',Y']} where {X'=2*X+rx}, {Y'=2*Y+ry}, and
    {rx,ry} range in {0..1}.
     
    The averaging will use the reliability weight of each sample
    {A[c,X',Y']}, assumed to {A[wch,X',Y']} if channel {wch} exists; if
    it does not, the reliability weight is assumed to be 1.0. If a
    sample {A[c,X',Y']} does not exist in {A}, its reliability weight is
    assumed to be zero. Any sample {A[c,X',Y']} with zero weight is
    considered completely unknown and thus excluded from the averaging.
    
    Apart from the size reduction, every output sample value {B[c,X,Y]}
    will be multiplied bt {scale}, for every channel {c} except {wch}.
    
    If {wch} is a valid channel index in {A}, the output weight sample
    {B[wch,X,Y]} will be a combination of the weights {A[wch,X',Y']}
    of the samples that are used to compute {B[c,X,Y]}.
    For this combination, {A[wch,X',Y']} is interpreted as the
    reciprocal of the estimated variance of the noise in the sample
    value {A[c,X',Y']}. 
    
    However, if any of those four pixel weights {A[wch,X',Y']} is zero,
    {B[wch,X,Y]} will be zero. This exception is justified by the
    assumption that a weight {A[wch,X',Y']=0} indicates that the
    quantity represented by {A[c,X',Y']} may not be integrable in the
    cell {[X',Y']} of {A}'s domain; in which case it will not be
    integrable in the cell {[X,Y]} of {B}'s domain.
    
    If the computed new weight {JG[2,X,Y]} is zero, or if any of the
    computed {JG[c,X,Y]} for any {c} is {NAN} or overflows, the value of
    {JG[0,X,Y]} and {JG[0,X,Y]} will be set to {NAN}, and {JG[2,X,Y]}
    will be set to zero. */

#endif
