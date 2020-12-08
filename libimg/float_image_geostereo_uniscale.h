#ifndef float_image_geostereo_uniscale_H
#define float_image_geostereo_uniscale_H

/* Tools for geometric stereo reconstruction from image pairs. */
/* Last edited on 2017-06-26 01:53:48 by stolfilocal */ 

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <float_image.h>

void float_image_geostereo_uniscale
  ( float_image_t *f1,   /* Image 1. */
    float_image_t *f2,   /* Image 2. */
    int32_t nwx,         /* Window width. */
    int32_t nwy,         /* Window height. */
    double dmin,         /* Minimum signed displacement (pixels). */
    double dmax,         /* Maximum signed displacement (pixels). */
    int32_t ncands,      /* Number of candidates to keep. */
    float_image_t **fdP, /* (OUT) Dispmap image. */
    float_image_t **fsP  /* (OUT) Scoremap image. */
  );
  /* Computes the disparity map for two images {f1,f2}.
  
    The two input images must have the same dimensions and channel
    counts, and must have been geometrically corrected so that rows with
    the same index are corresponding epipolar lines.
  
    The returned displacement map {fd} is an image with same dimensions
    as {img1} and {img2}, but with {ncands} channels. Each pixel of {fd}
    in column {x} and row {y} contains the best {ncands} displacement
    values found for the corresponding pixels of the two images.
    
    More precisley, a displacement value {d} for that pixel means that
    the local values around row {y} and columns {x+d} of {img1} seem to
    match those around row {y} and column {x-d} of {img2}.
    
    The procedure will try all displacements {d} with a granularity
    {1/3} pixel. The procedure considers only values {d} in the interval
    {[dmin _ dmax]}. The input images are interpolated when {d} is not
    an integer.
    
    The comparison is made using a sampling window with {nwx} columns
    and {nwy} rows centered at the points in question.  The two 
    dimensions must be odd.
    
    The procedure also returns a score map {fs}, with same dimensions
    and with {ncands} channels. Each sample of {fs} is the discrepancy
    {s} between the two image neighborhoods implied by the corresponding
    sample {d} from the displacement map. The discrepancy is defined as
    the average squared difference of the coresponding samples of {img1}
    and {img2} in the respective sampling windows.
    
    For each pixel of the image domain, the procedure enumerates all displacements
    of the form {i/3}, where {i} is an integer, in the range {[dmin _ dmax]};
    and computes their scores {sc[i]}.  Then it looks for local minima 
    in the list {sc[]}, and stores the {ncands} smallest minima in {fd,fs}.
    
    The images {fd,fs} are allocated by the procedure, and retuned in
    {*fdP} and {*fsP}, if these are not {NULL}. */

#endif
