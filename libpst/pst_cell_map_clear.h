#ifndef pst_cell_map_clear_H
#define pst_cell_map_clear_H

/* pst_cell_map_clear.h -- procedures for clearing regions of maps. */
/* Last edited on 2025-03-06 12:59:10 by stolfi */

#include <stdint.h>
#include <i2.h>
#include <vec.h>
#include <float_image.h>

void pst_cell_map_clear
  ( float_image_t *A,
    int32_t wch,
    i2_vec_t *imin,
    i2_vec_t *imax
  );
  /* The vectors {imin,imax} must have the same length {imin.ne=imax.ne}.
  
    For each {k} in {0..min.ne-1}, let {(xmin,ymin)=imin.e[k]} and {(xmax,ymax)=imax.e[k]}.
    The procedure sets to {NAN} every sample {A[c,x,y]} with {c!=wch}, {x} in {xmin..xmax},
    and {y} in {ymin..ymax}.  If {wch} is a valid channel index, also sets {A[wch,x,y]}
    to zero. */

#endif
