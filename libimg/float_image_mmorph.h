#ifndef float_image_mmorph_H
#define float_image_mmorph_H

/* Tools for mathematical morphology operations. */
/* Last edited on 2024-12-04 23:30:07 by stolfi */ 

#include <bool.h>
#include <r2.h>
#include <gauss_table.h>
#include <float_image.h>

float_image_t *float_image_mmorph_dilate(float_image_t *A, int32_t hw, double wt[]);
  /* Computes the dilation of image {A} by the tensor mask {wt'*wt} derived from
    the unidimensional weights {wt[0..2*hw]}.
  
    More precisely, the procedure computes an image {G} with the same
    size and channel count as {A}, such that each sample {G(c,x,y)} is
    
      {max { wt[dx+hw]*wt[dy+hw]*A[c,x+dx,y+dy] : dx,dy \in -hw..+hw }}
    
    where {A[c,x,y]} is {float_image_get_sample(A,c,x,y)}. The vector
    {wt} must be a unidimensional weight table with {2*hw+1} entries.
    Usually the entries are non-negative and symmetric around
    {wt[hw]}. */

/* !!! Add {float_image_mmorph_erode} !!! */

#endif
