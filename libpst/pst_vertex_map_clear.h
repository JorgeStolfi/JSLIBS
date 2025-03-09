#ifndef pst_vertex_map_clear_H
#define pst_vertex_map_clear_H

/* pst_vertex_map_clear.h -- procedures for clearing regions of maps. */
/* Last edited on 2025-03-06 12:58:47 by stolfi */

#include <stdint.h>
#include <i2.h>
#include <vec.h>
#include <float_image.h>

void pst_vertex_map_clear
  ( float_image_t *A,
    int32_t wch,
    i2_vec_t *imin,
    i2_vec_t *imax
  );
  /* The vectors {imin,imax} must have the same length {imin.ne=imax.ne}.
  
    For each {k} in {0..min.ne-1}, let {(xmin,ymin)=imin.e[k]} and {(xmax,ymax)=imax.e[k]}.
    The procedure sets to {NAN} every sample {A[c,x,y]} with {c!=wch}, {x} in {xmin..xmax+1},
    and {y} in {ymin..ymax+1}.  If {wch} is a valid channel index, also sets {A[wch,x,y]}
    to zero.
    
    Note that it assumes that the map elements are associated with the
    VERTICES of the integer grid, while the ranges {xmin..xmax} and
    {ymin..ymax} are assumed to be indices of CELLS. */

#endif
